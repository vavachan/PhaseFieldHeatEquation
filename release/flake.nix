{
  description = "OpenLB";

  inputs = {
    nixpkgs.url = github:NixOS/nixpkgs/nixos-24.11;
    nixpkgs-old.url = github:NixOS/nixpkgs/nixos-24.05;
  };

  outputs = { self, nixpkgs, nixpkgs-old, ... }: let
    system = "x86_64-linux";
    pkgs = import nixpkgs {
      inherit system;
      config.allowUnfree = true;
    };
    pkgs-old = import nixpkgs-old {
      inherit system;
      config.allowUnfree = true;
    };

    common-env = with pkgs;  [
    # common make dependencies
      gnumake

    # debugging
      gdb
      valgrind

    # formatting
      clang-tools

    # result presentation
      gnuplot
    ];

    mkShell = content: (pkgs.mkShell.override {
      stdenv = pkgs.stdenvNoCC;
    }) content;

  in {
    packages.${system} = {
      userguide = pkgs.stdenvNoCC.mkDerivation rec {
        name = "openlb-userguide";
        src = ./.;
        buildInputs = with pkgs; let
          custom-texlive = pkgs.texlive.combine {
            inherit (pkgs.texlive) scheme-small collection-langgerman latexmk xpatch xstring siunitx biblatex logreq palatino courier mathpazo helvetic multirow elsarticle widetable makecell pgfplots spath3 placeins abstract tocloft comment algorithmicx algorithms;
          };

        in [
          gnumake
          biber
          custom-texlive
        ];
        buildPhase = ''
          make userguide
        '';
        installPhase = ''
          mkdir -p $out
          cp doc/userGuide/olb-ug.pdf $out/
        '';
      };

      doxygen = pkgs.stdenvNoCC.mkDerivation rec {
        name = "openlb-doxygen";
        src = ./.;
        buildInputs = with pkgs; [
          gnumake
          doxygen
          graphviz
        ];
        buildPhase = ''
          make doxygen
        '';
        installPhase = ''
          mkdir -p $out
          cp -r doc/doxygen/html $out/
        '';
      };

      env-gcc = pkgs.mkShell {
        name = "openlb-env-gcc";
        buildInputs = common-env ++ (with pkgs; [
          gcc13
        ]);
        shellHook = ''
          export CXX=g++
          export CC=gcc

          export CXXFLAGS="-O3 -Wall -march=native -mtune=native -std=c++20"

          export PARALLEL_MODE=NONE

          export PLATFORMS="CPU_SISD"
        '';
      };

      env-emscripten = mkShell {
        name = "openlb-env-emscripten";
        buildInputs = common-env ++ (with pkgs; [
          gcc13
          emscripten
        ]);
        shellHook = ''
          export CXX=em++
          export CC=emcc

          export CXXFLAGS="-O3 -Wall -std=c++20 -Wno-missing-template-arg-list-after-template-kw"
          export LDFLAGS="-fPIC -g0 -s ALLOW_MEMORY_GROWTH -s MAXIMUM_MEMORY=4294967296 -s ENVIRONMENT=web -s EXPORTED_RUNTIME_METHODS=ccall,cwrap,UTF8ToString -fwasm-exceptions -lembind -s WASM=1"

          export PARALLEL_MODE=NONE
          export PLATFORMS="CPU_SISD"

          export FEATURES="EMSCRIPTEN"
        '';
      };

      env-gcc-openmpi = mkShell {
        name = "openlb-env-gcc-openmpi";
        buildInputs = common-env ++ (with pkgs; [
          gcc13
          mpi
        ]);
        shellHook = ''
          export CXX=mpic++
          export CC=gcc

          export CXXFLAGS="-O3 -g -Wall -march=native -mtune=native -std=c++20"

          export PARALLEL_MODE=MPI

          export PLATFORMS="CPU_SISD"
        '';
      };

      env-gcc-openmp = mkShell {
        name = "openlb-env-gcc-openmp";
        buildInputs = common-env ++ (with pkgs; [
          gcc13
        ]);
        shellHook = ''
          export CXX=g++
          export CC=gcc

          export CXXFLAGS="-O3 -Wall -march=native -mtune=native -std=c++20"

          export PARALLEL_MODE=OMP
          export OMPFLAGS="-fopenmp"

          export PLATFORMS="CPU_SISD"
        '';
      };

      env-gcc-openmpi-nlopt = mkShell {
        name = "openlb-env-nlopt";
        buildInputs = common-env ++ (with pkgs; [
          gcc13
          mpi
          vtk_9
          nlopt
        ]);
        shellHook = ''
          export CXX=mpic++
          export CC=gcc

          export CXXFLAGS="-O3 -Wall -march=native -mtune=native -std=c++20"

          export PARALLEL_MODE=MPI

          export PLATFORMS="CPU_SISD"

          export FEATURES="VTK"
          export VTK_VERSION=
          export CXXFLAGS="$CXXFLAGS -I${pkgs.vtk_9}/include/vtk -I${pkgs.nlopt}/include/nlopt"
          export LDFLAGS="$LDFLAGS -L${pkgs.vtk_9}/lib -L${pkgs.nlopt}/lib -lnlopt"
        '';
      };

      env-ibm-evaluation = pkgs.mkShell {
        name = "openlb-env-ibm-evaluation";
        buildInputs = common-env ++ (with pkgs; [
          gcc11
          gsl
          blas
        ]);
        shellHook = ''
          export CXX=g++
          export CC=gcc

          export CXXFLAGS="-O3 -Wall -march=native -mtune=native -std=c++20"
          export LDFLAGS="-lgsl -lblas"

          export PARALLEL_MODE=NONE

          export PLATFORMS="CPU_SISD"
        '';
      };

      env-clang = pkgs.mkShell.override {
        stdenv = pkgs.llvmPackages_19.libcxxStdenv;
      } {
        name = "openlb-env-clang";
        buildInputs = common-env ++ (with pkgs; [
          clang_19
          llvmPackages_19.bintools-unwrapped
          llvmPackages_19.openmp
        ]);
        shellHook = ''
          export CXX=clang++
          export CC=clang

          export CXXFLAGS="-O3 -Wall -march=native -mtune=native -std=c++20 -Wno-vla-cxx-extension"

          export PARALLEL_MODE=NONE

          export PLATFORMS="CPU_SISD"
        '';
      };

      env-cuda = mkShell {
        name = "openlb-env-cuda";
        buildInputs = common-env ++ (with pkgs; [
          gcc12
          cudaPackages.cudatoolkit
          cudaPackages.cuda_cudart
        ]);
        shellHook = ''
          export CXX=nvcc
          export CC=nvcc

          export PARALLEL_MODE=NONE

          export PLATFORMS="CPU_SISD GPU_CUDA"

          export CUDA_LDFLAGS="-L/run/opengl-driver/lib"
          # try to auto-fill CUDA_ARCH for first GPU (override when in doubt)
          export CUDA_ARCH=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n 1 | grep -o [0-9] | tr -d '\n')

          export FLOATING_POINT_TYPE=float

          export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/run/opengl-driver/lib
        '';
      };

      env-gcc-cuda = mkShell {
        name = "openlb-env-gcc-cuda";
        buildInputs = common-env ++ (with pkgs; [
          gcc12
          cudaPackages.cudatoolkit
          cudaPackages.cuda_cudart
        ]);
        shellHook = ''
          export CXX=g++
          export CC=gcc

          export CXXFLAGS="-O3 -Wall -march=native -mtune=native -std=c++20"

          export PARALLEL_MODE=NONE

          export PLATFORMS="CPU_SISD GPU_CUDA"

          export CUDA_CXX=nvcc
          export CUDA_CXXFLAGS="-O3 -std=c++20"
          export CUDA_LDFLAGS="-L/run/opengl-driver/lib"
          # try to auto-fill CUDA_ARCH for first GPU (override when in doubt)
          export CUDA_ARCH=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n 1 | grep -o [0-9] | tr -d '\n')

          export FLOATING_POINT_TYPE=float

          export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/run/opengl-driver/lib
        '';
      };

      env-gcc-openmpi-cuda = mkShell {
        name = "openlb-env-gcc-openmpi-cuda";
        buildInputs = common-env ++ (with pkgs; [
          gcc12
          cudaPackages.cudatoolkit
          cudaPackages.cuda_cudart
          (pkgs-old.mpi.override {
            cudaSupport = true;
          })
          vtk_9
        ]);
        shellHook = ''
          export CXX=mpic++
          export CC=gcc

          export CXXFLAGS="-O3 -g -Wall -march=native -mtune=native -std=c++20"

          export PARALLEL_MODE=MPI

          export PLATFORMS="CPU_SISD GPU_CUDA"

          export CUDA_CXX=nvcc
          export CUDA_CXXFLAGS="-O3 -std=c++20"
          export CUDA_LDFLAGS="-L/run/opengl-driver/lib"
          # try to auto-fill CUDA_ARCH for first GPU (override when in doubt)
          export CUDA_ARCH=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n 1 | grep -o [0-9] | tr -d '\n')

          export FLOATING_POINT_TYPE=float

          export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/run/opengl-driver/lib

          export FEATURES="VTK"
          export VTK_VERSION=
          export CXXFLAGS="$CXXFLAGS -I${pkgs.vtk_9}/include/vtk"
          export LDFLAGS="$LDFLAGS -L${pkgs.vtk_9}/lib"
        '';
      };

      env-heterogeneity-cpu = mkShell {
        name = "openlb-env-heterogeneity-cpu";
        buildInputs = common-env ++ (with pkgs; [
          gcc12
          mpi
        ]);
        shellHook = let
          arch = "tigerlake"; # match target CPU as closely as possible
        in ''
          export CXX=mpic++
          export CC=gcc

          export PLATFORMS="CPU_SISD CPU_SIMD"

          export CXXFLAGS="-O3 -Wall -march=${arch} -mtune=${arch} -std=c++20"

          export PARALLEL_MODE=HYBRID
          export OMPFLAGS="-fopenmp"
        '';
      };

      env-heterogeneity = mkShell {
        name = "openlb-env-heterogeneity";
        buildInputs = common-env ++ (with pkgs; [
          gcc12
          cudaPackages.cudatoolkit
          cudaPackages.cuda_cudart
          (pkgs-old.mpi.override {
            cudaSupport = true;
          })
        ]);
        shellHook = let
          arch = "tigerlake"; # match target CPU as closely as possible
        in ''
          export CXX=mpic++
          export CC=gcc

          export PLATFORMS="CPU_SISD CPU_SIMD GPU_CUDA"

          export CXXFLAGS="-O3 -Wall -march=${arch} -mtune=${arch} -std=c++20"

          export PARALLEL_MODE=HYBRID
          export OMPFLAGS="-fopenmp"

          export CUDA_CXX=nvcc
          export CUDA_CXXFLAGS="-O3 -std=c++20 --forward-unknown-to-host-compiler"
          export CUDA_LDFLAGS="-L/run/opengl-driver/lib"
          # try to auto-fill CUDA_ARCH for first GPU (override when in doubt)
          export CUDA_ARCH=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n 1 | grep -o [0-9] | tr -d '\n')

          export FLOATING_POINT_TYPE=float

          export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/run/opengl-driver/lib
        '';
      };

      env-fsi-cpu = mkShell {
        name = "openlb-env-fsi-cpu";
        buildInputs = common-env ++ (with pkgs; [
          gcc12
          mpi
          precice

          vtk_9

        # OpenVDB and its dependencies
          openvdb
          openexr
          boost
          tbb
        ]);
        shellHook = ''
          export CXX=mpic++
          export CC=gcc

          export CXXFLAGS="-O3 -Wall -march=native -mtune=native -std=c++20"

          export PARALLEL_MODE=HYBRID
          export OMPFLAGS="-fopenmp"

          export PLATFORMS="CPU_SISD"

          export FLOATING_POINT_TYPE=float

          export FEATURES="PRECICE VDB VTK"

          export VTK_VERSION=
          export CXXFLAGS="$CXXFLAGS -I${pkgs.vtk_9}/include/vtk"
          export LDFLAGS="$LDFLAGS -L${pkgs.vtk_9}/lib"
        '';
      };

      env-fsi = mkShell {
        name = "openlb-env-fsi";
        buildInputs = common-env ++ (with pkgs; [
          gcc12
          cudaPackages.cudatoolkit
          cudaPackages.cuda_cudart
          (pkgs-old.mpi.override {
            cudaSupport = true;
          })
          pkgs-old.precice

          vtk_9

        # OpenVDB and its dependencies
          openvdb
          openexr
          boost
          tbb
        ]);
        shellHook = ''
          export CXX=mpic++
          export CC=gcc

          export CXXFLAGS="-O3 -Wall -march=native -mtune=native -std=c++20"

          export PARALLEL_MODE=MPI

          export PLATFORMS="CPU_SISD GPU_CUDA"

          export CUDA_CXX=nvcc
          export CUDA_CXXFLAGS="-O3 -std=c++20"
          export CUDA_LDFLAGS="-L/run/opengl-driver/lib"
          # try to auto-fill CUDA_ARCH for first GPU (override when in doubt)
          export CUDA_ARCH=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n 1 | grep -o [0-9] | tr -d '\n')

          export FLOATING_POINT_TYPE=float

          export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/run/opengl-driver/lib

          export FEATURES="PRECICE VDB VTK"

          export VTK_VERSION=
          export CXXFLAGS="$CXXFLAGS -I${pkgs.vtk_9}/include/vtk"
          export LDFLAGS="$LDFLAGS -L${pkgs.vtk_9}/lib"
        '';
      };

      # Environment for code generation / CSE optimization
      env-code-generation = pkgs.mkShell {
        name = "openlb-env-code-generation";
        buildInputs = with pkgs; [
          gnumake
          gcc13
          gdb
          (python3.withPackages (python-pkgs: with python-pkgs; [
            mako
            sympy
          ]))
          (pkgs.writeShellScriptBin "generate-code" ''
            pushd script/codegen/cse
            make clean
            make TYPE=dynamics prepare_extraction
            make TYPE=dynamics extraction
            make TYPE=dynamics optimize_cse
            make TYPE=operator prepare_extraction
            make TYPE=operator extraction
            make TYPE=operator optimize_cse
            popd
          '')
        ];
        shellHook = ''
          export CXX=g++
          export CC=gcc

          export CXXFLAGS="-O1 -g -Wall -std=c++20"

          export PARALLEL_MODE=NONE

          export PLATFORMS="CPU_SISD"

          export FEATURES="INSPECT_DYNAMICS EXPORT_CODE_GENERATION_TARGETS"
        '';
      };
    };
  };
}
