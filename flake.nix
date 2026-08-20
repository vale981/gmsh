{
  description = "gmsh development environment";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs =
    { self, nixpkgs }:
    let
      forAllSystems = nixpkgs.lib.genAttrs [
        "x86_64-linux"
        "aarch64-linux"
        "x86_64-darwin"
        "aarch64-darwin"
      ];
    in
    {
      devShells = forAllSystems (
        system:
        let
          pkgs = nixpkgs.legacyPackages.${system};
          gmshDrv = pkgs.gmsh.overrideAttrs (_: {
            # Only used to inherit buildInputs/nativeBuildInputs below.
            doCheck = false;
          });
        in
        {
          default = pkgs.mkShell.override { stdenv = pkgs.llvmPackages.stdenv; } {
            inputsFrom = [ gmshDrv ];

            packages = with pkgs; [
              cmake
              ninja
              gfortran
              llvmPackages.clang-tools # clangd, clang-format, clang-tidy
              gdb
              ccache
              mmg
              (python3.withPackages (ps: [ ps.numpy ]))
            ];

            cmakeFlags = [
              "-DENABLE_BUILD_SHARED=ON"
              "-DENABLE_BUILD_DYNAMIC=ON"
              "-DENABLE_OPENMP=ON"
              "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON"
              "-DCMAKE_C_COMPILER_LAUNCHER=ccache"
              "-DCMAKE_CXX_COMPILER_LAUNCHER=ccache"
              "-DENABLE_MMG=ON"
              "-DMMG_INC=${pkgs.mmg.dev}/include"
              "-DMMG_LIB=${pkgs.mmg}/lib/libmmg.so"
              "-GNinja"
            ];

            shellHook = ''
              if [ -f build/compile_commands.json ] && [ ! -e compile_commands.json ]; then
                ln -s build/compile_commands.json compile_commands.json
              fi
              # gmsh.py locates libgmsh relative to its own directory (or via
              # ctypes.util.find_library, which ignores LD_LIBRARY_PATH on
              # Linux), so it never finds the lib in build/ on its own.
              if [ -f build/libgmsh.so.5.0 ] && [ ! -e api/libgmsh.so.5.0 ]; then
                ln -s ../build/libgmsh.so.5.0 api/libgmsh.so.5.0
              fi
              export PYTHONPATH="$PWD/api:$PYTHONPATH"
              export LD_LIBRARY_PATH="$PWD/build:$LD_LIBRARY_PATH"
            '';
          };
        }
      );
    };
}
