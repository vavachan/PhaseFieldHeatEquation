{
  description = "OpenLB";

  inputs = {
    nixpkgs.url = github:NixOS/nixpkgs/nixpkgs-unstable;
  };

  outputs = { self, nixpkgs, ... }: let
    system = "x86_64-linux";
    pkgs = import nixpkgs {
      inherit system;
      config.allowUnfree = true;
    };

    mkShell = content: (pkgs.mkShell.override {
      stdenv = pkgs.stdenvNoCC;
    }) content;

  in {
    packages.${system} = {
      # FHS for building virtualenv for code generation
      env-generate-code-fhs = let
        buildCommand = ''
          virtualenv tmp
          source tmp/bin/activate
          export CXXFLAGS="-std=c++20"
          export STDCXX=20
          export MAKE_NPROCS=$NIX_BUILD_CORES
          pip install --verbose "cppyy-cling==6.30.0" --no-binary=cppyy-cling --no-build-isolation
          pip install --verbose "cppyy-backend==1.15.2" --no-build-isolation
          pip install --verbose "CPyCppyy==1.12.16" --no-build-isolation
          pip install --verbose "cppyy==3.1.2" --no-build-isolation
          pip install "sympy==1.11.1"
          pip install "Mako==1.2.4"
        '';
      in pkgs.buildFHSUserEnv {
        name = "env-generate-code-fhs";
        targetPkgs = pkgs: (with pkgs; [
          python39
          python39Packages.pip
          python39Packages.virtualenv
          gnumake
          cmake
          zlib
          glibc
          glibc.dev
          gcc13
          (pkgs.writeScriptBin "build-virtualenv" buildCommand)
        ]);
        runScript = "bash";
      };

      # Bootstrap virtualenv containing cppyy for code generation (impure)
      env-generate-code-venv = pkgs.stdenvNoCC.mkDerivation rec {
        name = "openlb-env-generate-code-venv";
        buildInputs = with pkgs; [
          gnumake
          cmake
          zlib
          glibc
          glibc.dev
          gcc13
        ];
        buildCommand = ''
          originalVenvPath=$(pwd)/tmp
          mkdir $out
          ${self.packages.${system}.env-generate-code-fhs}/bin/env-generate-code-fhs "build-virtualenv"
          sed -i -e "s|$originalVenvPath|$out|g" tmp/bin/activate
          mv tmp/* $out/
        '';
      };

      # Call "make cse" in code generation environment
      generate-code = pkgs.buildFHSUserEnv {
        name = "openlb-generate-code";
        targetPkgs = pkgs: with pkgs; [
          gnumake
          cmake
          zlib
          glibc
          glibc.dev
          gcc13
          (pkgs.writeScriptBin "make-cse-in-codegen-env" ''
            #!/usr/bin/env bash
            source ${self.packages.${system}.env-generate-code-venv}/bin/activate
            IN_NIX_SHELL=1 make cse
          '')
        ];
        runScript = "make-cse-in-codegen-env";
      };
    };
  };
}
