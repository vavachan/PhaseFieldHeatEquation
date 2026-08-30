import numpy as np
import vtk
import xml.etree.ElementTree as ET
from pathlib import Path
from vtk.util.numpy_support import vtk_to_numpy



class VTK2Numpy:



    def __init__(self, path_pvd, other_input):
        self.path_pvd = path_pvd

        # Parse the PVD file
        tree = ET.parse(self.path_pvd)
        root = tree.getroot()

        # Dictionary to hold timestep and corresponding filename
        self.timestep_files = {}

        # Iterate through DataSet elements and extract timestep and file attributes
        for dataset in root.findall('.//DataSet'):
            timestep = dataset.attrib['timestep']
            file_name = dataset.attrib['file']
            self.timestep_files[timestep] = file_name

        if isinstance(other_input, str):
            self.global_metadata = self._process_global_metadata(other_input)

        elif isinstance(other_input, dict):
            if not set(other_input.keys()) == {'spacing', 'min', 'max', 'shape', 'base_points'}:
                raise ValueError('other_input is a dict, but has the wrong set of keys.')
            self.global_metadata = other_input

        else:
            raise TypeError("Argument other_input be either a str or a dict")



    def _get_vtk_multi_block(self, file_name):
        """Provides a handler to a VTK multi-block dataset.
        Args:
            filename: str: Path to the VTK file (extension .vtm)
        Returns:
            ...: vtk.vtkMultiBlockDataSet: Handler to the multi-block dataset
        """
        reader = vtk.vtkXMLMultiBlockDataReader()
        reader.SetFileName(file_name)
        reader.Update()

        return reader.GetOutput()



    def _read_vtk_data(self, top_level_multi_block, metadata, array_name):
        """Reads the content of a VTK multi-block dataset
        and returns to a dictionary (to make sure that the numbering is correct).
        Args:
            top_level_multi_block: vtk.vtkMultiBlockDataSet: Handler to a multi-block dataset
            array_name:            str:                      Name of the data array to extract
        Returns:
            {str, np.array}: Block data dictionary in the form {blockID, data}    
        """
        data_collected = {}
        for i in range(top_level_multi_block.GetNumberOfBlocks()):
            block = top_level_multi_block.GetBlock(i).GetBlock(0)
            data = vtk_to_numpy(block.GetPointData().GetArray(array_name))
            data_ext = data if len(data.shape)==3 else data.reshape((*data.shape, 1))

            data_list = []
            for j in range(metadata[i]['dimension']):
                data_list.append( data_ext[:,j].reshape(metadata[i]['shape']) )

            # shape: (dimension, z, y, x)
            data_collected[i] = np.stack(data_list)

        return data_collected



    def _read_vtk_metadata(self, top_level_multi_block, array_name):
        """Reads the metadata associated to a VTK multi-block dataset
        and returns to a dictionary (to make sure that the numbering is correct).
        The metadata entries consist of:
            dimension: int
            shape:     np.array(int)
            spacing:   np.float64
        Args:
            top_level_multi_block: vtk.vtkMultiBlockDataSet: Handler to a multi-block dataset
        Returns:
            {str, dict}: Block metadata dictionary in the form {blockID, metadata_for_that_block}    
        """
        metadata_collected = {}
        for i in range(top_level_multi_block.GetNumberOfBlocks()):
            block = top_level_multi_block.GetBlock(i).GetBlock(0)

            def get_dimension():
                return block.GetPointData().GetArray(array_name).GetNumberOfComponents()

            def get_spacing_and_check_uniqueness():
                spacings = set(np.float64(block.GetSpacing()))
                if not len(spacings) == 1:
                    raise ValueError(f"Error in the vtm file no. {i}. The spacings are not all the same.")
                return spacings.pop()

            def get_extent0_N():
                extent = block.GetExtent()
                return np.array([ extent[2*j] for j in range(3) ])[::-1]

            def get_extent1_N():
                extent = block.GetExtent()
                return np.array([ extent[2*j+1] for j in range(3) ])[::-1]

            def get_shape():
                return get_extent1_N() - get_extent0_N() + 1

            metadata_collected[i] = {
                'dimension':      get_dimension(),
                'shape':          get_shape(),
                'spacing':        get_spacing_and_check_uniqueness(),
            }

        return metadata_collected



    def _process_global_metadata (self, file_name):
        """Determines the metadata of the global computational domain.
        The metadata entries consist of:
            spacing:     np.float64
            min:         np.array(np.float64)
            max:         np.array(np.float64)
            shape:       np.array(int)
            base_points: {int, np.array(int)}
        Args:
            filename: str: Path to the VTK file (extension .vtm) representing the coordinates of the points of the computational domain
        returns:
            {str, dict}: Global metadata dictionary in the form {entry_name, entry_value}    
        """
        multi_block = self._get_vtk_multi_block(file_name)
        metadata = self._read_vtk_metadata(multi_block, 'coordinates')
        field_per_block = self._read_vtk_data(multi_block, metadata, 'coordinates')

        def get_spacing_and_check_uniqueness():
            spacings = set([ item['spacing'] for item in metadata.values() ])
            if not len(spacings) == 1:
                raise ValueError(f"Error in the vtm files. The spacings are not all the same.")
            return spacings.pop()

        def get_min():
            item = np.min( np.stack([ item[:, 0, 0, 0] for item in field_per_block.values() ]), axis=0 )[::-1]
            return item if item.shape[0]==3 else np.array((0, *item))

        def get_max():
            item = np.max( np.stack([ item[:,-1,-1,-1] for item in field_per_block.values() ]), axis=0 )[::-1]
            return item if item.shape[0]==3 else np.array((0, *item))

        def get_shape():
            return np.array(( get_max() - get_min() ) / get_spacing_and_check_uniqueness() + 1).astype(int)

        def get_base_points():
            base_points = {}

            for id_block, block in field_per_block.items():
                block_min = block[:,0,0,0 ][::-1]
                block_min_ext = block_min if block_min.shape[0]==3 else np.array((0, *block_min))
                base_points [id_block] = np.array( (block_min_ext - get_min()) / get_spacing_and_check_uniqueness() ).astype(int)

            return base_points

        return {
            'spacing':     get_spacing_and_check_uniqueness(),
            'min':         get_min(),
            'max':         get_max(),
            'shape':       get_shape(),
            'base_points': get_base_points(),
        }



    def process_field_by_index(self, index, field_name):
        return self.process_field_by_key (
            list(self.timestep_files.keys())[index],
            field_name
        )



    def process_field_by_key(self, key, field_name):
        return self.process_field_by_path (
            str(Path(self.path_pvd).parent) + '/' + self.timestep_files[key],
            field_name
        )



    def process_field_by_path(self, path_data, field_name):
        multi_block = self._get_vtk_multi_block(path_data)
        metadata = self._read_vtk_metadata(multi_block, field_name)
        field_per_block = self._read_vtk_data(multi_block, metadata, field_name)

        def get_dimension():
            dimensions = set([ item['dimension'] for item in metadata.values() ])
            if not len(dimensions) == 1:
                raise ValueError(f"Error in the vtm files. The dimensions are not all the same.")
            return dimensions.pop()
    
        global_data = np.zeros( np.array((get_dimension(), *self.global_metadata['shape'])), dtype=np.float64 )
    
        for i in field_per_block.keys():
            for j in range(get_dimension()):
                global_data[ j,
                    self.global_metadata['base_points'][i][0] : self.global_metadata['base_points'][i][0] + metadata[i]['shape'][0],
                    self.global_metadata['base_points'][i][1] : self.global_metadata['base_points'][i][1] + metadata[i]['shape'][1],
                    self.global_metadata['base_points'][i][2] : self.global_metadata['base_points'][i][2] + metadata[i]['shape'][2]
                ] = field_per_block[i][j,:,:,:]
        return global_data
