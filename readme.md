# Yocto/Pathtrace: Tiny Path Tracer

In this homework, I learned how to build a simple path tracer. In particular how to:

- write camera with depth of field,
- write a complex material,
- write a naive path tracer,
- write a path tracer with multiple importance sampling.

## Framework

The code uses the library [Yocto/GL](https://github.com/xelatihy/yocto-gl),
that is included in this project in the directory `yocto`.
We suggest to consult the documentation for the library that you can find
at the beginning of the header files. Also, since the library is getting improved
during the duration of the course, se suggest that you star it and watch it
on Github, so that you can notified as improvements are made.

In order to compile the code, you have to install
[Xcode](https://apps.apple.com/it/app/xcode/id497799835?mt=12)
on OsX, [Visual Studio 2019](https://visualstudio.microsoft.com/it/vs/) on Windows,
or a modern version of gcc or clang on Linux,
together with the tools [cmake](www.cmake.org) and [ninja](https://ninja-build.org).
The script `scripts/build.sh` will perform a simple build on OsX.
As discussed in class, we prefer to use
[Visual Studio Code](https://code.visualstudio.com), with
[C/C++](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools) and
[CMake Tools](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cmake-tools)
extensions, that we have configured to use for this course.

This repository also contains tests that are executed from the command line
as shown in `run.sh`. The rendered images are saved in the `out/` directory.
The results should match the ones in the directory `check/`.

## Base Functionality


- **Camera Sampling**
- **Naive Path tracing**
- **Brdf sampling**
- **Delta handling**
- **Light sampling**
- **Path tracing**

## Extra features

### Refraction

I implemented it by adding a new case in the switches of various functions for this new type of material, as suggested by the professor. Then I added the functions offered by yocto for refraction.

<img src="out/highres/18_refraction_1280_1024.jpg" alt="refraction image" style="width:720px;"/>

### Bilinear Patches

In this case I followed the code from [Cool Patches](https://link.springer.com/chapter/10.1007/978-1-4842-4427-2_8) for the intersection of the bilinear patch and then modified eval_normal and eval_position to account for this change. For the normals we can use the precalculated normals for each vertex and interpolate between them. We do the same for the positions but we interpolate using the vertices. For these interpolations we use the uv coordinates obtained from the patch intersection. These changes were made in the yocto_bvh, yocto_geometry and yocto_scene files.

### Denoising

As suggested by the professor I created a separate executable by compiling the code, and including OpenImageIO for handling jpeg and png files. The executable takes as input the image, the albedos and the normals and outputs the denoised image. At first I only used the color from material.color to generate the albedos, but following the [oidn documentation](https://www.openimagedenoise.org/documentation.html) I also tried to generate better ones by using the color behind the object, when handling transparent materials, and calculating the reflected color for delta reflective materials. This second approach gives better results, but we need more time to generate the inputs for the denoising.

#### Noisy

<img src="out/highres/11_bathroom1_1280_8.jpg" alt="noisy image" style="width:720px;"/>

#### Simple albedo

<img src="out/oidn-simple-albedo/11_bathroom1_720_256.jpg" alt="simple albedo denoised" style="width:720px;"/>

#### Better albedo

<img src="out/oidn-better-albedo/11_bathroom1_720_256.jpg" alt="better albedo denoised" style="width:720px;"/>

<div style="page-break-after: always;"></div>

### Hair Shading

For the hair rendering I followed pbrt's [implementation](https://www.pbrt.org/hair.pdf). I created separate files that contain almost all of the code for the hair material and integrated them in yocto by setting up the cmake files, all the files are in libs/yocto_extra. I also modified some code in yocto_scene and yocto_sceneio to better integrate the hair material and to read from the scene file the material properties.
For the BSDF calculations, since yocto uses world coordinates, I transform the normal, the outgoing and incoming directions to the coordinate system used by pbrt, which uses a frame that has as z direction the normal's direction and as y the tangent's direction. I also multiply for the cosine in the eval_hair_scattering function, following the yocto design.
The rendering test scene was created using a hair model found in [this repository](https://github.com/dsforza96/yocto-hair).

<img src="out/highres/19_hair_1280_1024.jpg" alt="hair rendering" style="width:720px;"/>
