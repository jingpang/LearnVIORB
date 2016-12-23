#ifndef MSGSYNCHRONIZER_H
#define MSGSYNCHRONIZER_H

#include <ros/ros.h>
#include <sensor_msgs/Image.h>
#include <queue>
#include <nav_msgs/Odometry.h>
#include <sensor_msgs/Imu.h>

using namespace std;

namespace ORBVIO
{
class MsgSynchronizer
{
public:
    enum Status{
        NOTINIT = 0,
        INIT,
        NORMAL
    };

    MsgSynchronizer(const double& imagedelay = 0.);
    ~MsgSynchronizer();

    // add messages in callbacks
    void addImageMsg(const sensor_msgs::ImageConstPtr &imgmsg);
    void addImuMsg(const sensor_msgs::ImuConstPtr &imumsg);

    // loop in main function to handle all messages
    bool getRecentMsgs(sensor_msgs::ImageConstPtr &imgmsg, std::vector<sensor_msgs::ImuConstPtr> &vimumsgs);

    void clearMsgs(void);

    // for message callback if needed
    void imageCallback(const sensor_msgs::ImageConstPtr& msg);
    void imuCallback(const sensor_msgs::ImuConstPtr& msg);

    //
    inline Status getStatus(void) {return _status;}

    double getImageDelaySec(void) const {return _imageMsgDelaySec;}

private:
    double _imageMsgDelaySec;  // image message delay to imu message, in seconds
    std::queue<sensor_msgs::ImageConstPtr> _imageMsgQueue;
    std::queue<sensor_msgs::ImuConstPtr> _imuMsgQueue;
    ros::Time _imuMsgTimeStart;
    Status _status;
    int _dataUnsyncCnt;
};

}

#endif // MSGSYNCHRONIZER_H
