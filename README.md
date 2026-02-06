How to run the application:

    Application is optimised and ran only on a Linux Machine. Make sure DGtal is installed on the Linux OS.

    Insert CMakeLists.txt into the main directory containing the build and Rice Grains folder (let's call it /app)

/app

    /app/RiceGrains
    /app/build
    /app/CMakeLists.txt

    Insert helloworld.cpp in /app/build directory

    Go to terminal and run the commands to build the application:

cd build cmake .. make

    Move the compiled application in /app/RiceGrains directory that contains all the images of rice grains

/app/RiceGrains

    /app/RiceGrains/Rice_japonais_seg_bin.pgm
    /app/RiceGrains/Rice_basmati_seg_bin.pgm
    /app/RiceGrains/Rice_camargue_seg_bin.pgm
    /app/RiceGrains/Rice_mixed2_seg_bin.pgm
    /app/RiceGrains/Rice_mixed3_seg_bin.pgm
    /app/RiceGrains/helloworld

    Run the next command to run the application:

cd .. cd RiceGrains ./helloworld
