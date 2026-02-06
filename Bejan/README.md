How to run the application:

0. Application is optimised and ran only on a Linux Machine. Make sure DGtal is installed on the Linux OS.

1. Insert CMakeLists.txt into the main directory containing the build and Rice Grains folder (let's call it /app)

/app
-	/app/RiceGrains
-	/app/build
-	/app/CMakeLists.txt

2. Insert helloworld.cpp in /app/build directory

3. Go to terminal and run the commands to build the application:

cd build
cmake ..
make

4. Move the compiled application in /app/RiceGrains directory that contains all the images of rice grains

/app/RiceGrains
-	/app/RiceGrains/Rice_japonais_seg_bin.pgm
-	/app/RiceGrains/Rice_basmati_seg_bin.pgm
-	/app/RiceGrains/Rice_camargue_seg_bin.pgm
-	/app/RiceGrains/Rice_mixed2_seg_bin.pgm
-	/app/RiceGrains/Rice_mixed3_seg_bin.pgm
-	/app/RiceGrains/helloworld

5. Run the next command to run the application:

cd ..
cd RiceGrains
./helloworld
