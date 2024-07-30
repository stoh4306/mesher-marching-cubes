- Build
  . mkdir build && cd build
  . cmake .. -G Ninja
  . ninja
- Test
  . unzip tsdf.7z
  . ./build/basic_marching_cubes tsdf.vox 0.5 mesh.obj
