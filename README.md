# RegistrationOptimization
Register two point clouds rigidly using closest point distances and levenberg-marquardt to optimize the transform

Reads in two point cloud vtk files and finds the best rigid transform that aligns the two clouds. The process is similar to ICP except a different set of transformation modes can be used e.g. nonrigid transforms since the optimization is based on a greedy hill climb algorithm.
