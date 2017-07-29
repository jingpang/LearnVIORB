FROM nvidia/cuda:8.0-devel-ubuntu14.04

# This Dockerfile attempts to run ORB_VIO
# Prerequisite:
# Currently, I could not find any easier method to manage GUI with Docker.
# 0. nvidia graphic card to handle GUI (X11 and OpenGL context) using nvidia's Docker image
# 1. Install docker and nvidia-docker
# https://github.com/NVIDIA/nvidia-docker
# 2. Modify ORB_VIO part in build.sh
# It would look like this:
# rosdep update
# source /opt/ros/$ROS_DISTRO/setup.bash
# export ROS_PACKAGE_PATH=${ROS_PACKAGE_PATH}:/home/${DOCKER_USER}/LearnVIORB/VIORB/Examples/ROS/ORB_VIO
# echo $ROS_PACKAGE_PATH
# cd Examples/ROS/ORB_VIO
# mkdir -p build
# cd build
# ..
# 3. Edit config/euroc.yaml
# E.g. test.InitVIOTmpPath: "/home/yourname/opensourcecode/OpenSourceORBVIO/tmp/"
# bagfile: "/home/yourname/LearnVIORB/dataset/V1_01_easy.bag"
# 4. You should have V1_01_easy.bag and mav0 under LearnVIORB directory
# 5. Run docker build . -t viorb
# 6. Then run ./run_docker.sh
# 7. I tried to make rosbag run automatically, but it still does not work.
# Here you need to run it manually. Open a new terminal and run
# docker ps (to get container id)
# docker exec -it <container_id> /bin/bash
# cd LearnVIORB
# ./run_rosbag

# set your username, uid, and gid
# you can run id command to get uid and gid
# from your host machine
# username also should be the same as your host machine
ENV DOCKER_USER=your_username
ENV UID=1002
ENV GID=1003

ENV ROS_DISTRO indigo
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8
ENV APP "glxgears"

RUN apt-key adv --keyserver hkp://p80.pool.sks-keyservers.net:80 --recv-keys 421C365BD9FF1F717815A3895523BAEEB01FA116
RUN echo "deb http://packages.ros.org/ros/ubuntu trusty main" > /etc/apt/sources.list.d/ros-latest.list
RUN apt-get -qq update && \
    apt-get -qq upgrade -y
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y \
    mesa-utils \
    google-perftools \
    libgoogle-perftools-dev \
    sudo \
    libglew-dev \
    cmake \
    libpython2.7-dev \
    libav-tools \
    libavcodec-dev \
    libavutil-dev \ 
    libavformat-dev \
    libswscale-dev \
    libjpeg-dev \ 
    libpng12-dev \
    libtiff5-dev \
    libopenexr-dev \
    git-core \
    libusb-1.0-0-dev \
    libblas-dev \
    liblapack-dev \
    qt5-default \
    libsuitesparse-dev \
    libboost-all-dev \
    build-essential \
    cmake \
    wget \
    unzip \
    git \
    libgtk2.0-dev \
    pkg-config \
    libavcodec-dev \
    libavformat-dev \
    libswscale-dev \
    python-dev \
    python-numpy \
    libtbb2 \
    libtbb-dev \
    libjpeg-dev \
    libpng-dev \
    libtiff-dev \
    libjasper-dev \
    libdc1394-22-dev \
    libcanberra-gtk-module \
    libdc1394-22 \
    python-rosdep \
    python-rosinstall \
    python-vcstools \
    ros-indigo-desktop-full

RUN rosdep init \
    && rosdep update

RUN groupadd --gid ${GID} ${DOCKER_USER}
RUN export uid=${UID} gid=${GID} && \
    mkdir -p /home/${DOCKER_USER} && \
    echo "${DOCKER_USER}:x:${uid}:${gid}:Teera,,,:/home/${DOCKER_USER}:/bin/bash" >> /etc/passwd && \
    echo "${DOCKER_USER}:x:${uid}:" >> /etc/group

RUN echo "${DOCKER_USER} ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/${DOCKER_USER} && \
    chmod 0440 /etc/sudoers.d/${DOCKER_USER} && \
    chown ${uid}:${gid} -R /home/${DOCKER_USER}

RUN wget https://sourceforge.net/projects/opencvlibrary/files/opencv-unix/2.4.13/opencv-2.4.13.zip && \
    unzip opencv-2.4.13.zip && \
    cd opencv-2.4.13/ && \
    mkdir build && \
    cd build && \
    cmake -D CMAKE_BUILD_TYPE=RELEASE -D CMAKE_INSTALL_PREFIX=/usr/local -D BUILD_PYTHON_SUPPORT=ON -D WITH_XINE=ON -D WITH_TBB=ON .. && \
    make -j$(nproc) && make install -j$(nproc) && \
    cd ../.. && \
    rm opencv-2.4.13.zip

RUN git clone https://github.com/ktossell/libuvc && \
    cd libuvc && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install

RUN git clone https://github.com/stevenlovegrove/Pangolin.git && \
    cd Pangolin && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make install

RUN wget http://bitbucket.org/eigen/eigen/get/3.2.10.zip && \
    unzip 3.2.10.zip && \
    cd eigen-eigen-b9cd8366d4e8 && \
    mkdir builddir && \
    cd builddir && \
    cmake .. && \
    make && \
    make install

RUN ln -s /usr/local/cuda/lib64/libcudart.so /usr/lib/libcudart.so
RUN ln -s /usr/local/cuda/lib64/libcudart.a /usr/lib/libcudart.a
RUN ln -s /usr/local/cuda/lib64/libcudart.so /usr/lib/libopencv_dep_cudart.so

ADD . /home/${DOCKER_USER}/LearnVIORB
ENV ROS_PACKAGE_PATH /home/${DOCKER_USER}/LearnVIORB/Examples/ROS/ORB_VIO
RUN ln -s /home/${DOCKER_USER}/LearnVIORB/Examples/ROS/ORB_VIO /opt/ros/${ROS_DISTRO}/share/ORB_VIO
RUN /bin/bash -c "source /opt/ros/$ROS_DISTRO/setup.bash ;\
                  cd /home/${DOCKER_USER}/LearnVIORB ;\
                  ./build.sh"

USER ${DOCKER_USER}
ENV HOME /home/${DOCKER_USER}
WORKDIR /home/${DOCKER_USER}
RUN sudo chown ${UID}:${GID} -R /home/${DOCKER_USER}

ADD run_rosbag /home/${DOCKER_USER}/run_rosbag
ADD run /home/${DOCKER_USER}/run
CMD /home/${DOCKER_USER}/run
