#  Vulkan Rendering

## Running the Project

* Install vulkan sdk
* Install glfw3

From the project root, run the following commands

1. Create a build directory

```mkdir build```

2. Go to build directory

```cd build```

3. Change the following paths

```VULKAN_ROOT_LOCATON="$HOME/Desktop/TUM/Courses/ss19/3d_scanning_and_motion_capture/" \
   VULKAN_SDK_VERSION="1.1.108.0" \
   VULKAN_SDK="$VULKAN_ROOT_LOCATON/vulkansdk-macos-$VULKAN_SDK_VERSION/macOS" \
   VK_ICD_FILENAMES="$VULKAN_SDK/etc/vulkan/icd.d/MoltenVK_icd.json" \
   VK_LAYER_PATH="$VULKAN_SDK/etc/vulkan/explicit_layer.d" \
   PATH="/usr/local/opt/python/libexec/bin:$VULKAN_SDK/bin:$PATH" \
   cmake ..
```

4. Make

```make```
 
5. Run

```./vulkan_rendering```