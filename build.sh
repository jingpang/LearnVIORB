echo "Configuring and building Thirdparty/DBoW2 ..."

cd Thirdparty/DBoW2
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j

cd ../../g2o

echo "Configuring and building Thirdparty/g2o ..."

mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j

cd ../../../

echo "Uncompress vocabulary ..."

#cd Vocabulary
#tar -xf ORBvoc.txt.tar.gz
#cd ..

echo "Configuring and building ORB_SLAM2 ..."

mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j
cd ..

echo "Build ROS node ..."

cd Examples/ROS/ORB_VIO
mkdir build
cd build
cmake .. -DROS_BUILD_TYPE=Release
make -j
cd ../../../../

echo ""
echo "Launch file in Examples/ROS/ORB_VIO/launch."
echo "Modify the configuration file config/euroc.yaml"
echo "Run as: roslaunch ORB_VIO testeuroc.launch"
echo ""

#echo "Converting vocabulary to binary"
#./tools/bin_vocabulary
