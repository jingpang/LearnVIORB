/**
* This file is part of ORB-SLAM2.
*
* Copyright (C) 2014-2016 Raúl Mur-Artal <raulmur at unizar dot es> (University of Zaragoza)
* For more information see <https://github.com/raulmur/ORB_SLAM2>
*
* ORB-SLAM2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* ORB-SLAM2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with ORB-SLAM2. If not, see <http://www.gnu.org/licenses/>.
*/

#include "LocalMapping.h"
#include "LoopClosing.h"
#include "ORBmatcher.h"
#include "Optimizer.h"

#include<mutex>
#include "Converter.h"

namespace ORB_SLAM2
{
using namespace std;

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
class KeyFrameInit
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    double mTimeStamp;
    KeyFrameInit* mpPrevKeyFrame;
    cv::Mat Twc;
    IMUPreintegrator mIMUPreInt;
    std::vector<IMUData> mvIMUData;
    Vector3d bg;


    KeyFrameInit(KeyFrame& kf):
        mTimeStamp(kf.mTimeStamp), mpPrevKeyFrame(NULL), Twc(kf.GetPoseInverse().clone()),
        mIMUPreInt(kf.GetIMUPreInt()), mvIMUData(kf.GetVectorIMUData()), bg(0,0,0)
    {
    }

    void ComputePreInt(void)
    {
        if(mpPrevKeyFrame == NULL)
        {
            return;
        }
        else
        {
            // Reset pre-integrator first
            mIMUPreInt.reset();

            if(mvIMUData.empty())
                return;

            // remember to consider the gap between the last KF and the first IMU
            {
                const IMUData& imu = mvIMUData.front();
                double dt = std::max(0., imu._t - mpPrevKeyFrame->mTimeStamp);
                mIMUPreInt.update(imu._g - bg,imu._a ,dt);  // Acc bias not considered here
            }
            // integrate each imu
            for(size_t i=0; i<mvIMUData.size(); i++)
            {
                const IMUData& imu = mvIMUData[i];
                double nextt;
                if(i==mvIMUData.size()-1)
                    nextt = mTimeStamp;         // last IMU, next is this KeyFrame
                else
                    nextt = mvIMUData[i+1]._t;  // regular condition, next is imu data

                // delta time
                double dt = std::max(0., nextt - imu._t);
                // update pre-integrator
                mIMUPreInt.update(imu._g - bg,imu._a ,dt);
            }
        }
    }

};

bool LocalMapping::GetUpdatingInitPoses(void)
{
    unique_lock<mutex> lock(mMutexUpdatingInitPoses);
    return mbUpdatingInitPoses;
}

void LocalMapping::SetUpdatingInitPoses(bool flag)
{
    unique_lock<mutex> lock(mMutexUpdatingInitPoses);
    mbUpdatingInitPoses = flag;
}

KeyFrame* LocalMapping::GetMapUpdateKF()
{
    unique_lock<mutex> lock(mMutexMapUpdateFlag);
    return mpMapUpdateKF;
}

bool LocalMapping::GetMapUpdateFlagForTracking()
{
    unique_lock<mutex> lock(mMutexMapUpdateFlag);
    return mbMapUpdateFlagForTracking;
}

void LocalMapping::SetMapUpdateFlagInTracking(bool bflag)
{
    unique_lock<mutex> lock(mMutexMapUpdateFlag);
    mbMapUpdateFlagForTracking = bflag;
    if(bflag)
    {
        mpMapUpdateKF = mpCurrentKeyFrame;
    }
}

bool LocalMapping::GetVINSInited(void)
{
    unique_lock<mutex> lock(mMutexVINSInitFlag);
    return mbVINSInited;
}

void LocalMapping::SetVINSInited(bool flag)
{
    unique_lock<mutex> lock(mMutexVINSInitFlag);
    mbVINSInited = flag;
}

bool LocalMapping::GetFirstVINSInited(void)
{
    unique_lock<mutex> lock(mMutexFirstVINSInitFlag);
    return mbFirstVINSInited;
}

void LocalMapping::SetFirstVINSInited(bool flag)
{
    unique_lock<mutex> lock(mMutexFirstVINSInitFlag);
    mbFirstVINSInited = flag;
}

cv::Mat LocalMapping::GetGravityVec()
{
    return mGravityVec.clone();
}

cv::Mat LocalMapping::GetRwiInit()
{
    return mRwiInit.clone();
}

void LocalMapping::VINSInitThread()
{
    unsigned long initedid = 0;
    cerr<<"start VINSInitThread"<<endl;
    while(1)
    {
        if(KeyFrame::nNextId > 2)
            if(!GetVINSInited() && mpCurrentKeyFrame->mnId > initedid)
            {
                initedid = mpCurrentKeyFrame->mnId;

                bool tmpbool = TryInitVIO();
                if(tmpbool)
                {
                    //SetFirstVINSInited(true);
                    //SetVINSInited(true);
                    break;
                }
            }
        usleep(3000);
        if(isFinished())
            break;
    }
}

bool LocalMapping::TryInitVIO(void)
{
    if(mpMap->KeyFramesInMap()<=mnLocalWindowSize)
        return false;

    static bool fopened = false;
    static ofstream fgw,fscale,fbiasa,fcondnum,ftime,fbiasg;
    string tmpfilepath = ConfigParam::getTmpFilePath();
    if(!fopened)
    {
        // Need to modify this to correct path
        fgw.open(tmpfilepath+"gw.txt");
        fscale.open(tmpfilepath+"scale.txt");
        fbiasa.open(tmpfilepath+"biasa.txt");
        fcondnum.open(tmpfilepath+"condnum.txt");
        ftime.open(tmpfilepath+"computetime.txt");
        fbiasg.open(tmpfilepath+"biasg.txt");
        if(fgw.is_open() && fscale.is_open() && fbiasa.is_open() &&
                fcondnum.is_open() && ftime.is_open() && fbiasg.is_open())
            fopened = true;
        else
        {
            cerr<<"file open error in TryInitVIO"<<endl;
            fopened = false;
        }
        fgw<<std::fixed<<std::setprecision(6);
        fscale<<std::fixed<<std::setprecision(6);
        fbiasa<<std::fixed<<std::setprecision(6);
        fcondnum<<std::fixed<<std::setprecision(6);
        ftime<<std::fixed<<std::setprecision(6);
        fbiasg<<std::fixed<<std::setprecision(6);
    }

    Optimizer::GlobalBundleAdjustemnt(mpMap, 10);

    // Extrinsics
    cv::Mat Tbc = ConfigParam::GetMatTbc();
    cv::Mat Rbc = Tbc.rowRange(0,3).colRange(0,3);
    cv::Mat pbc = Tbc.rowRange(0,3).col(3);
    cv::Mat Rcb = Rbc.t();
    cv::Mat pcb = -Rcb*pbc;

    if(ConfigParam::GetRealTimeFlag())
    {
        // Wait KeyFrame Culling.
        // 1. if KeyFrame Culling is running, wait until finished.
        // 2. if KFs are being copied, then don't run KeyFrame Culling (in KeyFrameCulling function)
        while(GetFlagCopyInitKFs())
        {
            usleep(1000);
        }
    }
    SetFlagCopyInitKFs(true);

    // Use all KeyFrames in map to compute
    vector<KeyFrame*> vScaleGravityKF = mpMap->GetAllKeyFrames();
    int N = vScaleGravityKF.size();
    KeyFrame* pNewestKF = vScaleGravityKF[N-1];
    vector<cv::Mat> vTwc;
    vector<IMUPreintegrator> vIMUPreInt;
    // Store initialization-required KeyFrame data
    vector<KeyFrameInit*> vKFInit;

    for(int i=0;i<N;i++)
    {
        KeyFrame* pKF = vScaleGravityKF[i];
        vTwc.push_back(pKF->GetPoseInverse());
        vIMUPreInt.push_back(pKF->GetIMUPreInt());
        KeyFrameInit* pkfi = new KeyFrameInit (*pKF);
        if(i!=0)
        {
            pkfi->mpPrevKeyFrame = vKFInit[i-1];
        }
        vKFInit.push_back(pkfi);
    }

    SetFlagCopyInitKFs(false);

    // Step 1.
    // Try to compute initial gyro bias, using optimization with Gauss-Newton
    Vector3d bgest = Optimizer::OptimizeInitialGyroBias(vTwc,vIMUPreInt);
    //Vector3d bgest = Optimizer::OptimizeInitialGyroBias(vScaleGravityKF);

    // Update biasg and pre-integration in LocalWindow. Remember to reset back to zero
    for(int i=0;i<N;i++)
    {
        vKFInit[i]->bg = bgest;
    }
    for(int i=0;i<N;i++)
    {
        vKFInit[i]->ComputePreInt();
    }

    // Solve A*x=B for x=[s,gw] 4x1 vector
    cv::Mat A = cv::Mat::zeros(3*(N-2),4,CV_32F);
    cv::Mat B = cv::Mat::zeros(3*(N-2),1,CV_32F);
    cv::Mat I3 = cv::Mat::eye(3,3,CV_32F);

    // Step 2.
    // Approx Scale and Gravity vector in 'world' frame (first KF's camera frame)
    for(int i=0; i<N-2; i++)
    {
        //KeyFrameInit* pKF1 = vKFInit[i];//vScaleGravityKF[i];
        KeyFrameInit* pKF2 = vKFInit[i+1];
        KeyFrameInit* pKF3 = vKFInit[i+2];
        // Delta time between frames
        double dt12 = pKF2->mIMUPreInt.getDeltaTime();
        double dt23 = pKF3->mIMUPreInt.getDeltaTime();
        // Pre-integrated measurements
        cv::Mat dp12 = Converter::toCvMat(pKF2->mIMUPreInt.getDeltaP());
        cv::Mat dv12 = Converter::toCvMat(pKF2->mIMUPreInt.getDeltaV());
        cv::Mat dp23 = Converter::toCvMat(pKF3->mIMUPreInt.getDeltaP());

        // Pose of camera in world frame
        cv::Mat Twc1 = vTwc[i].clone();//pKF1->GetPoseInverse();
        cv::Mat Twc2 = vTwc[i+1].clone();//pKF2->GetPoseInverse();
        cv::Mat Twc3 = vTwc[i+2].clone();//pKF3->GetPoseInverse();
        // Position of camera center
        cv::Mat pc1 = Twc1.rowRange(0,3).col(3);
        cv::Mat pc2 = Twc2.rowRange(0,3).col(3);
        cv::Mat pc3 = Twc3.rowRange(0,3).col(3);
        // Rotation of camera, Rwc
        cv::Mat Rc1 = Twc1.rowRange(0,3).colRange(0,3);
        cv::Mat Rc2 = Twc2.rowRange(0,3).colRange(0,3);
        cv::Mat Rc3 = Twc3.rowRange(0,3).colRange(0,3);

        // Stack to A/B matrix
        // lambda*s + beta*g = gamma
        cv::Mat lambda = (pc2-pc1)*dt23 + (pc2-pc3)*dt12;
        cv::Mat beta = 0.5*I3*(dt12*dt12*dt23 + dt12*dt23*dt23);
        cv::Mat gamma = (Rc3-Rc2)*pcb*dt12 + (Rc1-Rc2)*pcb*dt23 + Rc1*Rcb*dp12*dt23 - Rc2*Rcb*dp23*dt12 - Rc1*Rcb*dv12*dt12*dt23;
        lambda.copyTo(A.rowRange(3*i+0,3*i+3).col(0));
        beta.copyTo(A.rowRange(3*i+0,3*i+3).colRange(1,4));
        gamma.copyTo(B.rowRange(3*i+0,3*i+3));
        // Tested the formulation in paper, -gamma. Then the scale and gravity vector is -xx

        // Debug log
        //cout<<"iter "<<i<<endl;
    }
    // Use svd to compute A*x=B, x=[s,gw] 4x1 vector
    // A = u*w*vt, u*w*vt*x=B
    // Then x = vt'*winv*u'*B
    cv::Mat w,u,vt;
    // Note w is 4x1 vector by SVDecomp()
    // A is changed in SVDecomp() with cv::SVD::MODIFY_A for speed
    cv::SVDecomp(A,w,u,vt,cv::SVD::MODIFY_A);
    // Debug log
    //cout<<"u:"<<endl<<u<<endl;
    //cout<<"vt:"<<endl<<vt<<endl;
    //cout<<"w:"<<endl<<w<<endl;

    // Compute winv
    cv::Mat winv=cv::Mat::eye(4,4,CV_32F);
    for(int i=0;i<4;i++)
    {
        if(fabs(w.at<float>(i))<1e-10)
        {
            w.at<float>(i) += 1e-10;
            // Test log
            cerr<<"w(i) < 1e-10, w="<<endl<<w<<endl;
        }

        winv.at<float>(i,i) = 1./w.at<float>(i);
    }
    // Then x = vt'*winv*u'*B
    cv::Mat x = vt.t()*winv*u.t()*B;

    // x=[s,gw] 4x1 vector
    double sstar = x.at<float>(0);    // scale should be positive
    cv::Mat gwstar = x.rowRange(1,4);   // gravity should be about ~9.8

    // Debug log
    //cout<<"scale sstar: "<<sstar<<endl;
    //cout<<"gwstar: "<<gwstar.t()<<", |gwstar|="<<cv::norm(gwstar)<<endl;

    // Test log
    if(w.type()!=I3.type() || u.type()!=I3.type() || vt.type()!=I3.type())
        cerr<<"different mat type, I3,w,u,vt: "<<I3.type()<<","<<w.type()<<","<<u.type()<<","<<vt.type()<<endl;

    // Step 3.
    // Use gravity magnitude 9.8 as constraint
    // gI = [0;0;1], the normalized gravity vector in an inertial frame, NED type with no orientation.
    cv::Mat gI = cv::Mat::zeros(3,1,CV_32F);
    gI.at<float>(2) = 1;
    // Normalized approx. gravity vecotr in world frame
    cv::Mat gwn = gwstar/cv::norm(gwstar);
    // Debug log
    //cout<<"gw normalized: "<<gwn<<endl;

    // vhat = (gI x gw) / |gI x gw|
    cv::Mat gIxgwn = gI.cross(gwn);
    double normgIxgwn = cv::norm(gIxgwn);
    cv::Mat vhat = gIxgwn/normgIxgwn;
    double theta = std::atan2(normgIxgwn,gI.dot(gwn));
    // Debug log
    //cout<<"vhat: "<<vhat<<", theta: "<<theta*180.0/M_PI<<endl;

    Eigen::Vector3d vhateig = Converter::toVector3d(vhat);
    Eigen::Matrix3d RWIeig = Sophus::SO3::exp(vhateig*theta).matrix();
    cv::Mat Rwi = Converter::toCvMat(RWIeig);
    cv::Mat GI = gI*ConfigParam::GetG();//9.8012;
    // Solve C*x=D for x=[s,dthetaxy,ba] (1+2+3)x1 vector
    cv::Mat C = cv::Mat::zeros(3*(N-2),6,CV_32F);
    cv::Mat D = cv::Mat::zeros(3*(N-2),1,CV_32F);

    for(int i=0; i<N-2; i++)
    {
        //KeyFrameInit* pKF1 = vKFInit[i];//vScaleGravityKF[i];
        KeyFrameInit* pKF2 = vKFInit[i+1];
        KeyFrameInit* pKF3 = vKFInit[i+2];
        // Delta time between frames
        double dt12 = pKF2->mIMUPreInt.getDeltaTime();
        double dt23 = pKF3->mIMUPreInt.getDeltaTime();
        // Pre-integrated measurements
        cv::Mat dp12 = Converter::toCvMat(pKF2->mIMUPreInt.getDeltaP());
        cv::Mat dv12 = Converter::toCvMat(pKF2->mIMUPreInt.getDeltaV());
        cv::Mat dp23 = Converter::toCvMat(pKF3->mIMUPreInt.getDeltaP());
        cv::Mat Jpba12 = Converter::toCvMat(pKF2->mIMUPreInt.getJPBiasa());
        cv::Mat Jvba12 = Converter::toCvMat(pKF2->mIMUPreInt.getJVBiasa());
        cv::Mat Jpba23 = Converter::toCvMat(pKF3->mIMUPreInt.getJPBiasa());
        // Pose of camera in world frame
        cv::Mat Twc1 = vTwc[i].clone();//pKF1->GetPoseInverse();
        cv::Mat Twc2 = vTwc[i+1].clone();//pKF2->GetPoseInverse();
        cv::Mat Twc3 = vTwc[i+2].clone();//pKF3->GetPoseInverse();
        // Position of camera center
        cv::Mat pc1 = Twc1.rowRange(0,3).col(3);
        cv::Mat pc2 = Twc2.rowRange(0,3).col(3);
        cv::Mat pc3 = Twc3.rowRange(0,3).col(3);
        // Rotation of camera, Rwc
        cv::Mat Rc1 = Twc1.rowRange(0,3).colRange(0,3);
        cv::Mat Rc2 = Twc2.rowRange(0,3).colRange(0,3);
        cv::Mat Rc3 = Twc3.rowRange(0,3).colRange(0,3);
        // Stack to C/D matrix
        // lambda*s + phi*dthetaxy + zeta*ba = psi
        cv::Mat lambda = (pc2-pc1)*dt23 + (pc2-pc3)*dt12;
        cv::Mat phi = - 0.5*(dt12*dt12*dt23 + dt12*dt23*dt23)*Rwi*SkewSymmetricMatrix(GI);  // note: this has a '-', different to paper
        cv::Mat zeta = Rc2*Rcb*Jpba23*dt12 + Rc1*Rcb*Jvba12*dt12*dt23 - Rc1*Rcb*Jpba12*dt23;
        cv::Mat psi = (Rc1-Rc2)*pcb*dt23 + Rc1*Rcb*dp12*dt23 - (Rc2-Rc3)*pcb*dt12
                     - Rc2*Rcb*dp23*dt12 - Rc1*Rcb*dv12*dt23*dt12 - 0.5*Rwi*GI*(dt12*dt12*dt23 + dt12*dt23*dt23); // note:  - paper
        lambda.copyTo(C.rowRange(3*i+0,3*i+3).col(0));
        phi.colRange(0,2).copyTo(C.rowRange(3*i+0,3*i+3).colRange(1,3)); //only the first 2 columns, third term in dtheta is zero, here compute dthetaxy 2x1.
        zeta.copyTo(C.rowRange(3*i+0,3*i+3).colRange(3,6));
        psi.copyTo(D.rowRange(3*i+0,3*i+3));

        // Debug log
        //cout<<"iter "<<i<<endl;
    }

    // Use svd to compute C*x=D, x=[s,dthetaxy,ba] 6x1 vector
    // C = u*w*vt, u*w*vt*x=D
    // Then x = vt'*winv*u'*D
    cv::Mat w2,u2,vt2;
    // Note w2 is 6x1 vector by SVDecomp()
    // C is changed in SVDecomp() with cv::SVD::MODIFY_A for speed
    cv::SVDecomp(C,w2,u2,vt2,cv::SVD::MODIFY_A);
    // Debug log
    //cout<<"u2:"<<endl<<u2<<endl;
    //cout<<"vt2:"<<endl<<vt2<<endl;
    //cout<<"w2:"<<endl<<w2<<endl;

    // Compute winv
    cv::Mat w2inv=cv::Mat::eye(6,6,CV_32F);
    for(int i=0;i<6;i++)
    {
        if(fabs(w2.at<float>(i))<1e-10)
        {
            w2.at<float>(i) += 1e-10;
            // Test log
            cerr<<"w2(i) < 1e-10, w="<<endl<<w2<<endl;
        }

        w2inv.at<float>(i,i) = 1./w2.at<float>(i);
    }
    // Then y = vt'*winv*u'*D
    cv::Mat y = vt2.t()*w2inv*u2.t()*D;

    double s_ = y.at<float>(0);
    cv::Mat dthetaxy = y.rowRange(1,3);
    cv::Mat dbiasa_ = y.rowRange(3,6);
    Vector3d dbiasa_eig = Converter::toVector3d(dbiasa_);

    // dtheta = [dx;dy;0]
    cv::Mat dtheta = cv::Mat::zeros(3,1,CV_32F);
    dthetaxy.copyTo(dtheta.rowRange(0,2));
    Eigen::Vector3d dthetaeig = Converter::toVector3d(dtheta);
    // Rwi_ = Rwi*exp(dtheta)
    Eigen::Matrix3d Rwieig_ = RWIeig*Sophus::SO3::exp(dthetaeig).matrix();
    cv::Mat Rwi_ = Converter::toCvMat(Rwieig_);


    // Debug log
    {
        cv::Mat gwbefore = Rwi*GI;
        cv::Mat gwafter = Rwi_*GI;
        cout<<"Time: "<<mpCurrentKeyFrame->mTimeStamp - mnStartTime<<", sstar: "<<sstar<<", s: "<<s_<<endl;

        fgw<<mpCurrentKeyFrame->mTimeStamp<<" "
           <<gwafter.at<float>(0)<<" "<<gwafter.at<float>(1)<<" "<<gwafter.at<float>(2)<<" "
           <<gwbefore.at<float>(0)<<" "<<gwbefore.at<float>(1)<<" "<<gwbefore.at<float>(2)<<" "
           <<endl;
        fscale<<mpCurrentKeyFrame->mTimeStamp<<" "
              <<s_<<" "<<sstar<<" "<<endl;
        fbiasa<<mpCurrentKeyFrame->mTimeStamp<<" "
              <<dbiasa_.at<float>(0)<<" "<<dbiasa_.at<float>(1)<<" "<<dbiasa_.at<float>(2)<<" "<<endl;
        fcondnum<<mpCurrentKeyFrame->mTimeStamp<<" "
                <<w2.at<float>(0)<<" "<<w2.at<float>(1)<<" "<<w2.at<float>(2)<<" "<<w2.at<float>(3)<<" "
                <<w2.at<float>(4)<<" "<<w2.at<float>(5)<<" "<<endl;
        //        ftime<<mpCurrentKeyFrame->mTimeStamp<<" "
        //             <<(t3-t0)/cv::getTickFrequency()*1000<<" "<<endl;
        fbiasg<<mpCurrentKeyFrame->mTimeStamp<<" "
              <<bgest(0)<<" "<<bgest(1)<<" "<<bgest(2)<<" "<<endl;

        ofstream fRwi(tmpfilepath+"Rwi.txt");
        fRwi<<Rwieig_(0,0)<<" "<<Rwieig_(0,1)<<" "<<Rwieig_(0,2)<<" "
            <<Rwieig_(1,0)<<" "<<Rwieig_(1,1)<<" "<<Rwieig_(1,2)<<" "
            <<Rwieig_(2,0)<<" "<<Rwieig_(2,1)<<" "<<Rwieig_(2,2)<<endl;
        fRwi.close();
    }


    // ********************************
    // Todo:
    // Add some logic or strategy to confirm init status
    bool bVIOInited = false;
    if(mbFirstTry)
    {
        mbFirstTry = false;
        mnStartTime = mpCurrentKeyFrame->mTimeStamp;
    }
    if(pNewestKF->mTimeStamp - mnStartTime >= ConfigParam::GetVINSInitTime())
    {
        bVIOInited = true;
    }

    if(bVIOInited)
    {
        // Set NavState , scale and bias for all KeyFrames
        // Scale
        double scale = s_;
        mnVINSInitScale = s_;
        // gravity vector in world frame
        cv::Mat gw = Rwi_*GI;
        mGravityVec = gw.clone();
        Vector3d gweig = Converter::toVector3d(gw);
        mRwiInit = Rwi_.clone();

        // Update NavState for the KeyFrames not in vScaleGravityKF
        // Update Tcw-type pose for these KeyFrames, need mutex lock
        if(ConfigParam::GetRealTimeFlag())
        {
            // Stop local mapping
            RequestStop();

            // Wait until Local Mapping has effectively stopped
            while(!isStopped() && !isFinished())
            {
                usleep(1000);
            }
        }

        SetUpdatingInitPoses(true);
        {
            unique_lock<mutex> lock(mpMap->mMutexMapUpdate);

            int cnt=0;
            for(vector<KeyFrame*>::const_iterator vit=vScaleGravityKF.begin(), vend=vScaleGravityKF.end(); vit!=vend; vit++,cnt++)
            {
                KeyFrame* pKF = *vit;
                if(pKF->isBad()) continue;
                if(pKF!=vScaleGravityKF[cnt]) cerr<<"pKF!=vScaleGravityKF[cnt], id: "<<pKF->mnId<<" != "<<vScaleGravityKF[cnt]->mnId<<endl;
                // Position and rotation of visual SLAM
                cv::Mat wPc = pKF->GetPoseInverse().rowRange(0,3).col(3);                   // wPc
                cv::Mat Rwc = pKF->GetPoseInverse().rowRange(0,3).colRange(0,3);            // Rwc
                // Set position and rotation of navstate
                cv::Mat wPb = scale*wPc + Rwc*pcb;
                pKF->SetNavStatePos(Converter::toVector3d(wPb));
                pKF->SetNavStateRot(Converter::toMatrix3d(Rwc*Rcb));
                // Update bias of Gyr & Acc
                pKF->SetNavStateBiasGyr(bgest);
                pKF->SetNavStateBiasAcc(dbiasa_eig);
                // Set delta_bias to zero. (only updated during optimization)
                pKF->SetNavStateDeltaBg(Eigen::Vector3d::Zero());
                pKF->SetNavStateDeltaBa(Eigen::Vector3d::Zero());
                // Step 4.
                // compute velocity
                if(pKF != vScaleGravityKF.back())
                {
                    KeyFrame* pKFnext = pKF->GetNextKeyFrame();
                    if(!pKFnext) cerr<<"pKFnext is NULL, cnt="<<cnt<<", pKFnext:"<<pKFnext<<endl;
                    if(pKFnext!=vScaleGravityKF[cnt+1]) cerr<<"pKFnext!=vScaleGravityKF[cnt+1], cnt="<<cnt<<", id: "<<pKFnext->mnId<<" != "<<vScaleGravityKF[cnt+1]->mnId<<endl;
                    // IMU pre-int between pKF ~ pKFnext
                    const IMUPreintegrator& imupreint = pKFnext->GetIMUPreInt();
                    // Time from this(pKF) to next(pKFnext)
                    double dt = imupreint.getDeltaTime();                                       // deltaTime
                    cv::Mat dp = Converter::toCvMat(imupreint.getDeltaP());       // deltaP
                    cv::Mat Jpba = Converter::toCvMat(imupreint.getJPBiasa());    // J_deltaP_biasa
                    cv::Mat wPcnext = pKFnext->GetPoseInverse().rowRange(0,3).col(3);           // wPc next
                    cv::Mat Rwcnext = pKFnext->GetPoseInverse().rowRange(0,3).colRange(0,3);    // Rwc next

                    cv::Mat vel = - 1./dt*( scale*(wPc - wPcnext) + (Rwc - Rwcnext)*pcb + Rwc*Rcb*(dp + Jpba*dbiasa_) + 0.5*gw*dt*dt );
                    Eigen::Vector3d veleig = Converter::toVector3d(vel);
                    pKF->SetNavStateVel(veleig);
                }
                else
                {
                    cerr<<"-----------here is the last KF in vScaleGravityKF------------"<<endl;
                    // If this is the last KeyFrame, no 'next' KeyFrame exists
                    KeyFrame* pKFprev = pKF->GetPrevKeyFrame();
                    if(!pKFprev) cerr<<"pKFprev is NULL, cnt="<<cnt<<endl;
                    if(pKFprev!=vScaleGravityKF[cnt-1]) cerr<<"pKFprev!=vScaleGravityKF[cnt-1], cnt="<<cnt<<", id: "<<pKFprev->mnId<<" != "<<vScaleGravityKF[cnt-1]->mnId<<endl;
                    const IMUPreintegrator& imupreint_prev_cur = pKF->GetIMUPreInt();
                    double dt = imupreint_prev_cur.getDeltaTime();
                    Eigen::Matrix3d Jvba = imupreint_prev_cur.getJVBiasa();
                    Eigen::Vector3d dv = imupreint_prev_cur.getDeltaV();
                    //
                    Eigen::Vector3d velpre = pKFprev->GetNavState().Get_V();
                    Eigen::Matrix3d rotpre = pKFprev->GetNavState().Get_RotMatrix();
                    Eigen::Vector3d veleig = velpre + gweig*dt + rotpre*( dv + Jvba*dbiasa_eig );
                    pKF->SetNavStateVel(veleig);
                }
            }

            // Re-compute IMU pre-integration at last. Should after usage of pre-int measurements.
            for(vector<KeyFrame*>::const_iterator vit=vScaleGravityKF.begin(), vend=vScaleGravityKF.end(); vit!=vend; vit++)
            {
                KeyFrame* pKF = *vit;
                if(pKF->isBad()) continue;
                pKF->ComputePreInt();
            }

            // Update poses (multiply metric scale)
            vector<KeyFrame*> mspKeyFrames = mpMap->GetAllKeyFrames();
            for(std::vector<KeyFrame*>::iterator sit=mspKeyFrames.begin(), send=mspKeyFrames.end(); sit!=send; sit++)
            {
                KeyFrame* pKF = *sit;
                cv::Mat Tcw = pKF->GetPose();
                cv::Mat tcw = Tcw.rowRange(0,3).col(3)*scale;
                tcw.copyTo(Tcw.rowRange(0,3).col(3));
                pKF->SetPose(Tcw);
            }
            vector<MapPoint*> mspMapPoints = mpMap->GetAllMapPoints();
            for(std::vector<MapPoint*>::iterator sit=mspMapPoints.begin(), send=mspMapPoints.end(); sit!=send; sit++)
            {
                MapPoint* pMP = *sit;
                //pMP->SetWorldPos(pMP->GetWorldPos()*scale);
                pMP->UpdateScale(scale);
            }
            std::cout<<std::endl<<"... Map scale updated ..."<<std::endl<<std::endl;

            // Update NavStates
            if(pNewestKF!=mpCurrentKeyFrame)
            {
                KeyFrame* pKF;

                // step1. bias&d_bias
                pKF = pNewestKF;
                do
                {
                    pKF = pKF->GetNextKeyFrame();

                    // Update bias of Gyr & Acc
                    pKF->SetNavStateBiasGyr(bgest);
                    pKF->SetNavStateBiasAcc(dbiasa_eig);
                    // Set delta_bias to zero. (only updated during optimization)
                    pKF->SetNavStateDeltaBg(Eigen::Vector3d::Zero());
                    pKF->SetNavStateDeltaBa(Eigen::Vector3d::Zero());
                }while(pKF!=mpCurrentKeyFrame);

                // step2. re-compute pre-integration
                pKF = pNewestKF;
                do
                {
                    pKF = pKF->GetNextKeyFrame();

                    pKF->ComputePreInt();
                }while(pKF!=mpCurrentKeyFrame);

                // step3. update pos/rot
                pKF = pNewestKF;
                do
                {
                    pKF = pKF->GetNextKeyFrame();

                    // Update rot/pos
                    // Position and rotation of visual SLAM
                    cv::Mat wPc = pKF->GetPoseInverse().rowRange(0,3).col(3);                   // wPc
                    cv::Mat Rwc = pKF->GetPoseInverse().rowRange(0,3).colRange(0,3);            // Rwc
                    cv::Mat wPb = wPc + Rwc*pcb;
                    pKF->SetNavStatePos(Converter::toVector3d(wPb));
                    pKF->SetNavStateRot(Converter::toMatrix3d(Rwc*Rcb));

                    //pKF->SetNavState();

                    if(pKF != mpCurrentKeyFrame)
                    {
                        KeyFrame* pKFnext = pKF->GetNextKeyFrame();
                        // IMU pre-int between pKF ~ pKFnext
                        const IMUPreintegrator& imupreint = pKFnext->GetIMUPreInt();
                        // Time from this(pKF) to next(pKFnext)
                        double dt = imupreint.getDeltaTime();                                       // deltaTime
                        cv::Mat dp = Converter::toCvMat(imupreint.getDeltaP());       // deltaP
                        cv::Mat Jpba = Converter::toCvMat(imupreint.getJPBiasa());    // J_deltaP_biasa
                        cv::Mat wPcnext = pKFnext->GetPoseInverse().rowRange(0,3).col(3);           // wPc next
                        cv::Mat Rwcnext = pKFnext->GetPoseInverse().rowRange(0,3).colRange(0,3);    // Rwc next

                        cv::Mat vel = - 1./dt*( (wPc - wPcnext) + (Rwc - Rwcnext)*pcb + Rwc*Rcb*(dp + Jpba*dbiasa_) + 0.5*gw*dt*dt );
                        Eigen::Vector3d veleig = Converter::toVector3d(vel);
                        pKF->SetNavStateVel(veleig);
                    }
                    else
                    {
                        // If this is the last KeyFrame, no 'next' KeyFrame exists
                        KeyFrame* pKFprev = pKF->GetPrevKeyFrame();
                        const IMUPreintegrator& imupreint_prev_cur = pKF->GetIMUPreInt();
                        double dt = imupreint_prev_cur.getDeltaTime();
                        Eigen::Matrix3d Jvba = imupreint_prev_cur.getJVBiasa();
                        Eigen::Vector3d dv = imupreint_prev_cur.getDeltaV();
                        //
                        Eigen::Vector3d velpre = pKFprev->GetNavState().Get_V();
                        Eigen::Matrix3d rotpre = pKFprev->GetNavState().Get_RotMatrix();
                        Eigen::Vector3d veleig = velpre + gweig*dt + rotpre*( dv + Jvba*dbiasa_eig );
                        pKF->SetNavStateVel(veleig);
                    }

                }while(pKF!=mpCurrentKeyFrame);

            }

            std::cout<<std::endl<<"... Map NavState updated ..."<<std::endl<<std::endl;

            SetFirstVINSInited(true);
            SetVINSInited(true);
        }
        SetUpdatingInitPoses(false);

        if(ConfigParam::GetRealTimeFlag())
        {
            Release();
        }

        // Run global BA after inited
        unsigned long nGBAKF = mpCurrentKeyFrame->mnId;
        //Optimizer::GlobalBundleAdjustmentNavState(mpMap,mGravityVec,10,NULL,nGBAKF,false);
        Optimizer::GlobalBundleAdjustmentNavStatePRV(mpMap,mGravityVec,10,NULL,nGBAKF,false);
        cerr<<"finish global BA after vins init"<<endl;

        if(ConfigParam::GetRealTimeFlag())
        {
            // Update pose
            // Stop local mapping, and
            RequestStop();

            // Wait until Local Mapping has effectively stopped
            while(!isStopped() && !isFinished())
            {
                usleep(1000);
            }


            cv::Mat cvTbc = ConfigParam::GetMatTbc();

            {
                unique_lock<mutex> lock(mpMap->mMutexMapUpdate);

                // Correct keyframes starting at map first keyframe
                list<KeyFrame*> lpKFtoCheck(mpMap->mvpKeyFrameOrigins.begin(),mpMap->mvpKeyFrameOrigins.end());

                while(!lpKFtoCheck.empty())
                {
                    KeyFrame* pKF = lpKFtoCheck.front();
                    const set<KeyFrame*> sChilds = pKF->GetChilds();
                    cv::Mat Twc = pKF->GetPoseInverse();
                    for(set<KeyFrame*>::const_iterator sit=sChilds.begin();sit!=sChilds.end();sit++)
                    {
                        KeyFrame* pChild = *sit;
                        if(pChild->mnBAGlobalForKF!=nGBAKF)
                        {
                            cerr<<"correct KF after gBA in VI init: "<<pChild->mnId<<endl;
                            cv::Mat Tchildc = pChild->GetPose()*Twc;
                            pChild->mTcwGBA = Tchildc*pKF->mTcwGBA;//*Tcorc*pKF->mTcwGBA;
                            pChild->mnBAGlobalForKF=nGBAKF;

                            // Set NavStateGBA and correct the P/V/R
                            pChild->mNavStateGBA = pChild->GetNavState();
                            cv::Mat TwbGBA = Converter::toCvMatInverse(cvTbc*pChild->mTcwGBA);
                            Matrix3d RwbGBA = Converter::toMatrix3d(TwbGBA.rowRange(0,3).colRange(0,3));
                            Vector3d PwbGBA = Converter::toVector3d(TwbGBA.rowRange(0,3).col(3));
                            Matrix3d Rw1 = pChild->mNavStateGBA.Get_RotMatrix();
                            Vector3d Vw1 = pChild->mNavStateGBA.Get_V();
                            Vector3d Vw2 = RwbGBA*Rw1.transpose()*Vw1;   // bV1 = bV2 ==> Rwb1^T*wV1 = Rwb2^T*wV2 ==> wV2 = Rwb2*Rwb1^T*wV1
                            pChild->mNavStateGBA.Set_Pos(PwbGBA);
                            pChild->mNavStateGBA.Set_Rot(RwbGBA);
                            pChild->mNavStateGBA.Set_Vel(Vw2);
                        }
                        lpKFtoCheck.push_back(pChild);
                    }

                    pKF->mTcwBefGBA = pKF->GetPose();
                    //pKF->SetPose(pKF->mTcwGBA);
                    pKF->mNavStateBefGBA = pKF->GetNavState();
                    pKF->SetNavState(pKF->mNavStateGBA);
                    pKF->UpdatePoseFromNS(cvTbc);

                    lpKFtoCheck.pop_front();

                }

                // Correct MapPoints
                const vector<MapPoint*> vpMPs = mpMap->GetAllMapPoints();

                for(size_t i=0; i<vpMPs.size(); i++)
                {
                    MapPoint* pMP = vpMPs[i];

                    if(pMP->isBad())
                        continue;

                    if(pMP->mnBAGlobalForKF==nGBAKF)
                    {
                        // If optimized by Global BA, just update
                        pMP->SetWorldPos(pMP->mPosGBA);
                    }
                    else
                    {
                        // Update according to the correction of its reference keyframe
                        KeyFrame* pRefKF = pMP->GetReferenceKeyFrame();

                        if(pRefKF->mnBAGlobalForKF!=nGBAKF)
                            continue;

                        // Map to non-corrected camera
                        cv::Mat Rcw = pRefKF->mTcwBefGBA.rowRange(0,3).colRange(0,3);
                        cv::Mat tcw = pRefKF->mTcwBefGBA.rowRange(0,3).col(3);
                        cv::Mat Xc = Rcw*pMP->GetWorldPos()+tcw;

                        // Backproject using corrected camera
                        cv::Mat Twc = pRefKF->GetPoseInverse();
                        cv::Mat Rwc = Twc.rowRange(0,3).colRange(0,3);
                        cv::Mat twc = Twc.rowRange(0,3).col(3);

                        pMP->SetWorldPos(Rwc*Xc+twc);
                    }
                }

                cout << "Map updated!" << endl;

                // Map updated, set flag for Tracking
                SetMapUpdateFlagInTracking(true);

                // Release LocalMapping

                Release();
            }
        }

        SetFlagInitGBAFinish(true);
    }

    for(int i=0;i<N;i++)
    {
        if(vKFInit[i])
            delete vKFInit[i];
    }

    return bVIOInited;
}

void LocalMapping::AddToLocalWindow(KeyFrame* pKF)
{
    mlLocalKeyFrames.push_back(pKF);
    if(mlLocalKeyFrames.size() > mnLocalWindowSize)
    {
        mlLocalKeyFrames.pop_front();
    }
    else
    {
        KeyFrame* pKF0 = mlLocalKeyFrames.front();
        while(mlLocalKeyFrames.size() < mnLocalWindowSize && pKF0->GetPrevKeyFrame()!=NULL)
        {
            pKF0 = pKF0->GetPrevKeyFrame();
            mlLocalKeyFrames.push_front(pKF0);
        }
    }
}

void LocalMapping::DeleteBadInLocalWindow(void)
{
    std::list<KeyFrame*>::iterator lit = mlLocalKeyFrames.begin();
    while(lit != mlLocalKeyFrames.end())
    {
        KeyFrame* pKF = *lit;
        //Test log
        if(!pKF) cout<<"pKF null?"<<endl;
        if(pKF->isBad())
        {
            lit = mlLocalKeyFrames.erase(lit);
        }
        else
        {
            lit++;
        }
    }
}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

LocalMapping::LocalMapping(Map *pMap, const float bMonocular, ConfigParam* pParams):
    mbMonocular(bMonocular), mbResetRequested(false), mbFinishRequested(false), mbFinished(true), mpMap(pMap),
    mbAbortBA(false), mbStopped(false), mbStopRequested(false), mbNotStop(false), mbAcceptKeyFrames(true)
{
    mpParams = pParams;
    mnLocalWindowSize = ConfigParam::GetLocalWindowSize();
    cout<<"mnLocalWindowSize:"<<mnLocalWindowSize<<endl;

    mbVINSInited = false;
    mbFirstTry = true;
    mbFirstVINSInited = false;

    mbUpdatingInitPoses = false;
    mbCopyInitKFs = false;
    mbInitGBAFinish = false;
}

void LocalMapping::SetLoopCloser(LoopClosing* pLoopCloser)
{
    mpLoopCloser = pLoopCloser;
}

void LocalMapping::SetTracker(Tracking *pTracker)
{
    mpTracker=pTracker;
}

void LocalMapping::Run()
{

    mbFinished = false;

    while(1)
    {
        // Tracking will see that Local Mapping is busy
        SetAcceptKeyFrames(false);

        // Check if there are keyframes in the queue
        if(CheckNewKeyFrames())
        {
            // Local Window also updated in below function
            // BoW conversion and insertion in Map
            ProcessNewKeyFrame();

            // Check recent MapPoints
            MapPointCulling();

            // Triangulate new MapPoints
            CreateNewMapPoints();

            if(!CheckNewKeyFrames())
            {
                // Find more matches in neighbor keyframes and fuse point duplications
                SearchInNeighbors();
            }

            mbAbortBA = false;

            if(!CheckNewKeyFrames() && !stopRequested())
            {
                // Local BA
                if(mpMap->KeyFramesInMap()>2)
                {
                    if(!GetVINSInited())
                    {
                        //Optimizer::LocalBundleAdjustment(mpCurrentKeyFrame,mlLocalKeyFrames,&mbAbortBA, mpMap, this);
                        Optimizer::LocalBundleAdjustment(mpCurrentKeyFrame,&mbAbortBA,mpMap,this);
                    }
                    else
                    {
                        //Optimizer::LocalBundleAdjustmentNavStatePRV(mpCurrentKeyFrame,mlLocalKeyFrames,&mbAbortBA, mpMap, mGravityVec, this);
                        Optimizer::LocalBAPRVIDP(mpCurrentKeyFrame,mlLocalKeyFrames,&mbAbortBA, mpMap, mGravityVec, this);
                    }
                }

                // Visual-Inertial initialization for non-realtime mode
                if(!ConfigParam::GetRealTimeFlag())
                {
                    // Try to initialize VIO, if not inited
                    if(!GetVINSInited())
                    {
                        bool tmpbool = TryInitVIO();
                        SetVINSInited(tmpbool);
                        if(tmpbool)
                        {
                            // Update map scale
                            mpMap->UpdateScale(mnVINSInitScale);
                            // Set initialization flag
                            SetFirstVINSInited(true);
                        }
                    }
                }

                // May set bad for KF in LocalWindow
                // Check redundant local Keyframes
                KeyFrameCulling();
            }

            if(GetFlagInitGBAFinish())
                mpLoopCloser->InsertKeyFrame(mpCurrentKeyFrame);
        }
        else if(Stop())
        {
            // Safe area to stop
            while(isStopped() && !CheckFinish())
            {
                usleep(3000);
            }
            if(CheckFinish())
                break;
        }

        ResetIfRequested();

        // Tracking will see that Local Mapping is busy
        SetAcceptKeyFrames(true);

        if(CheckFinish())
            break;

        usleep(3000);
    }

    SetFinish();
}

void LocalMapping::InsertKeyFrame(KeyFrame *pKF)
{
    unique_lock<mutex> lock(mMutexNewKFs);
    mlNewKeyFrames.push_back(pKF);
    mbAbortBA=true;
}


bool LocalMapping::CheckNewKeyFrames()
{
    unique_lock<mutex> lock(mMutexNewKFs);
    return(!mlNewKeyFrames.empty());
}

void LocalMapping::ProcessNewKeyFrame()
{
    {
        unique_lock<mutex> lock(mMutexNewKFs);
        mpCurrentKeyFrame = mlNewKeyFrames.front();
        mlNewKeyFrames.pop_front();
    }

    // Compute Bags of Words structures
    mpCurrentKeyFrame->ComputeBoW();

    // Associate MapPoints to the new keyframe and update normal and descriptor
    const vector<MapPoint*> vpMapPointMatches = mpCurrentKeyFrame->GetMapPointMatches();

    for(size_t i=0; i<vpMapPointMatches.size(); i++)
    {
        MapPoint* pMP = vpMapPointMatches[i];
        if(pMP)
        {
            if(!pMP->isBad())
            {
                if(!pMP->IsInKeyFrame(mpCurrentKeyFrame))
                {
                    pMP->AddObservation(mpCurrentKeyFrame, i);
                    pMP->UpdateNormalAndDepth();
                    pMP->ComputeDistinctiveDescriptors();
                }
                else // this can only happen for new stereo points inserted by the Tracking
                {
                    mlpRecentAddedMapPoints.push_back(pMP);
                }
            }
        }
    }    

    // Update links in the Covisibility Graph
    mpCurrentKeyFrame->UpdateConnections();

    // Delete bad KF in LocalWindow
    DeleteBadInLocalWindow();
    // Add Keyframe to LocalWindow
    AddToLocalWindow(mpCurrentKeyFrame);

    // Insert Keyframe in Map
    mpMap->AddKeyFrame(mpCurrentKeyFrame);
}

void LocalMapping::MapPointCulling()
{
    // Check Recent Added MapPoints
    list<MapPoint*>::iterator lit = mlpRecentAddedMapPoints.begin();
    const unsigned long int nCurrentKFid = mpCurrentKeyFrame->mnId;

    int nThObs;
    if(mbMonocular)
        nThObs = 2;
    else
        nThObs = 3;
    const int cnThObs = nThObs;

    while(lit!=mlpRecentAddedMapPoints.end())
    {
        MapPoint* pMP = *lit;
        if(pMP->isBad())
        {
            lit = mlpRecentAddedMapPoints.erase(lit);
        }
        else if(pMP->GetFoundRatio()<0.25f )
        {
            pMP->SetBadFlag();
            lit = mlpRecentAddedMapPoints.erase(lit);
        }
        else if(((int)nCurrentKFid-(int)pMP->mnFirstKFid)>=2 && pMP->Observations()<=cnThObs)
        {
            pMP->SetBadFlag();
            lit = mlpRecentAddedMapPoints.erase(lit);
        }
        else if(((int)nCurrentKFid-(int)pMP->mnFirstKFid)>=3)
            lit = mlpRecentAddedMapPoints.erase(lit);
        else
            lit++;
    }
}

void LocalMapping::CreateNewMapPoints()
{
    // Retrieve neighbor keyframes in covisibility graph
    int nn = 10;
    if(mbMonocular)
        nn=20;
    const vector<KeyFrame*> vpNeighKFs = mpCurrentKeyFrame->GetBestCovisibilityKeyFrames(nn);

    ORBmatcher matcher(0.6,false);

    cv::Mat Rcw1 = mpCurrentKeyFrame->GetRotation();
    cv::Mat Rwc1 = Rcw1.t();
    cv::Mat tcw1 = mpCurrentKeyFrame->GetTranslation();
    cv::Mat Tcw1(3,4,CV_32F);
    Rcw1.copyTo(Tcw1.colRange(0,3));
    tcw1.copyTo(Tcw1.col(3));
    cv::Mat Ow1 = mpCurrentKeyFrame->GetCameraCenter();

    const float &fx1 = mpCurrentKeyFrame->fx;
    const float &fy1 = mpCurrentKeyFrame->fy;
    const float &cx1 = mpCurrentKeyFrame->cx;
    const float &cy1 = mpCurrentKeyFrame->cy;
    const float &invfx1 = mpCurrentKeyFrame->invfx;
    const float &invfy1 = mpCurrentKeyFrame->invfy;

    const float ratioFactor = 1.5f*mpCurrentKeyFrame->mfScaleFactor;

    int nnew=0;

    // Search matches with epipolar restriction and triangulate
    for(size_t i=0; i<vpNeighKFs.size(); i++)
    {
        if(i>0 && CheckNewKeyFrames())
            return;

        KeyFrame* pKF2 = vpNeighKFs[i];

        // Check first that baseline is not too short
        cv::Mat Ow2 = pKF2->GetCameraCenter();
        cv::Mat vBaseline = Ow2-Ow1;
        const float baseline = cv::norm(vBaseline);

        if(!mbMonocular)
        {
            if(baseline<pKF2->mb)
            continue;
        }
        else
        {
            const float medianDepthKF2 = pKF2->ComputeSceneMedianDepth(2);
            const float ratioBaselineDepth = baseline/medianDepthKF2;

            if(ratioBaselineDepth<0.01)
                continue;
        }

        // Compute Fundamental Matrix
        cv::Mat F12 = ComputeF12(mpCurrentKeyFrame,pKF2);

        // Search matches that fullfil epipolar constraint
        vector<pair<size_t,size_t> > vMatchedIndices;
        matcher.SearchForTriangulation(mpCurrentKeyFrame,pKF2,F12,vMatchedIndices,false);

        cv::Mat Rcw2 = pKF2->GetRotation();
        cv::Mat Rwc2 = Rcw2.t();
        cv::Mat tcw2 = pKF2->GetTranslation();
        cv::Mat Tcw2(3,4,CV_32F);
        Rcw2.copyTo(Tcw2.colRange(0,3));
        tcw2.copyTo(Tcw2.col(3));

        const float &fx2 = pKF2->fx;
        const float &fy2 = pKF2->fy;
        const float &cx2 = pKF2->cx;
        const float &cy2 = pKF2->cy;
        const float &invfx2 = pKF2->invfx;
        const float &invfy2 = pKF2->invfy;

        // Triangulate each match
        const int nmatches = vMatchedIndices.size();
        for(int ikp=0; ikp<nmatches; ikp++)
        {
            const int &idx1 = vMatchedIndices[ikp].first;
            const int &idx2 = vMatchedIndices[ikp].second;

            const cv::KeyPoint &kp1 = mpCurrentKeyFrame->mvKeysUn[idx1];
            const float kp1_ur=mpCurrentKeyFrame->mvuRight[idx1];
            bool bStereo1 = kp1_ur>=0;

            const cv::KeyPoint &kp2 = pKF2->mvKeysUn[idx2];
            const float kp2_ur = pKF2->mvuRight[idx2];
            bool bStereo2 = kp2_ur>=0;

            // Check parallax between rays
            cv::Mat xn1 = (cv::Mat_<float>(3,1) << (kp1.pt.x-cx1)*invfx1, (kp1.pt.y-cy1)*invfy1, 1.0);
            cv::Mat xn2 = (cv::Mat_<float>(3,1) << (kp2.pt.x-cx2)*invfx2, (kp2.pt.y-cy2)*invfy2, 1.0);

            cv::Mat ray1 = Rwc1*xn1;
            cv::Mat ray2 = Rwc2*xn2;
            const float cosParallaxRays = ray1.dot(ray2)/(cv::norm(ray1)*cv::norm(ray2));

            float cosParallaxStereo = cosParallaxRays+1;
            float cosParallaxStereo1 = cosParallaxStereo;
            float cosParallaxStereo2 = cosParallaxStereo;

            if(bStereo1)
                cosParallaxStereo1 = cos(2*atan2(mpCurrentKeyFrame->mb/2,mpCurrentKeyFrame->mvDepth[idx1]));
            else if(bStereo2)
                cosParallaxStereo2 = cos(2*atan2(pKF2->mb/2,pKF2->mvDepth[idx2]));

            cosParallaxStereo = min(cosParallaxStereo1,cosParallaxStereo2);

            cv::Mat x3D;
            if(cosParallaxRays<cosParallaxStereo && cosParallaxRays>0 && (bStereo1 || bStereo2 || cosParallaxRays<0.9998))
            {
                // Linear Triangulation Method
                cv::Mat A(4,4,CV_32F);
                A.row(0) = xn1.at<float>(0)*Tcw1.row(2)-Tcw1.row(0);
                A.row(1) = xn1.at<float>(1)*Tcw1.row(2)-Tcw1.row(1);
                A.row(2) = xn2.at<float>(0)*Tcw2.row(2)-Tcw2.row(0);
                A.row(3) = xn2.at<float>(1)*Tcw2.row(2)-Tcw2.row(1);

                cv::Mat w,u,vt;
                cv::SVD::compute(A,w,u,vt,cv::SVD::MODIFY_A| cv::SVD::FULL_UV);

                x3D = vt.row(3).t();

                if(x3D.at<float>(3)==0)
                    continue;

                // Euclidean coordinates
                x3D = x3D.rowRange(0,3)/x3D.at<float>(3);

            }
            else if(bStereo1 && cosParallaxStereo1<cosParallaxStereo2)
            {
                x3D = mpCurrentKeyFrame->UnprojectStereo(idx1);                
            }
            else if(bStereo2 && cosParallaxStereo2<cosParallaxStereo1)
            {
                x3D = pKF2->UnprojectStereo(idx2);
            }
            else
                continue; //No stereo and very low parallax

            cv::Mat x3Dt = x3D.t();

            //Check triangulation in front of cameras
            float z1 = Rcw1.row(2).dot(x3Dt)+tcw1.at<float>(2);
            if(z1<=0)
                continue;

            float z2 = Rcw2.row(2).dot(x3Dt)+tcw2.at<float>(2);
            if(z2<=0)
                continue;

            //Check reprojection error in first keyframe
            const float &sigmaSquare1 = mpCurrentKeyFrame->mvLevelSigma2[kp1.octave];
            const float x1 = Rcw1.row(0).dot(x3Dt)+tcw1.at<float>(0);
            const float y1 = Rcw1.row(1).dot(x3Dt)+tcw1.at<float>(1);
            const float invz1 = 1.0/z1;

            if(!bStereo1)
            {
                float u1 = fx1*x1*invz1+cx1;
                float v1 = fy1*y1*invz1+cy1;
                float errX1 = u1 - kp1.pt.x;
                float errY1 = v1 - kp1.pt.y;
                if((errX1*errX1+errY1*errY1)>5.991*sigmaSquare1)
                    continue;
            }
            else
            {
                float u1 = fx1*x1*invz1+cx1;
                float u1_r = u1 - mpCurrentKeyFrame->mbf*invz1;
                float v1 = fy1*y1*invz1+cy1;
                float errX1 = u1 - kp1.pt.x;
                float errY1 = v1 - kp1.pt.y;
                float errX1_r = u1_r - kp1_ur;
                if((errX1*errX1+errY1*errY1+errX1_r*errX1_r)>7.8*sigmaSquare1)
                    continue;
            }

            //Check reprojection error in second keyframe
            const float sigmaSquare2 = pKF2->mvLevelSigma2[kp2.octave];
            const float x2 = Rcw2.row(0).dot(x3Dt)+tcw2.at<float>(0);
            const float y2 = Rcw2.row(1).dot(x3Dt)+tcw2.at<float>(1);
            const float invz2 = 1.0/z2;
            if(!bStereo2)
            {
                float u2 = fx2*x2*invz2+cx2;
                float v2 = fy2*y2*invz2+cy2;
                float errX2 = u2 - kp2.pt.x;
                float errY2 = v2 - kp2.pt.y;
                if((errX2*errX2+errY2*errY2)>5.991*sigmaSquare2)
                    continue;
            }
            else
            {
                float u2 = fx2*x2*invz2+cx2;
                float u2_r = u2 - mpCurrentKeyFrame->mbf*invz2;
                float v2 = fy2*y2*invz2+cy2;
                float errX2 = u2 - kp2.pt.x;
                float errY2 = v2 - kp2.pt.y;
                float errX2_r = u2_r - kp2_ur;
                if((errX2*errX2+errY2*errY2+errX2_r*errX2_r)>7.8*sigmaSquare2)
                    continue;
            }

            //Check scale consistency
            cv::Mat normal1 = x3D-Ow1;
            float dist1 = cv::norm(normal1);

            cv::Mat normal2 = x3D-Ow2;
            float dist2 = cv::norm(normal2);

            if(dist1==0 || dist2==0)
                continue;

            const float ratioDist = dist2/dist1;
            const float ratioOctave = mpCurrentKeyFrame->mvScaleFactors[kp1.octave]/pKF2->mvScaleFactors[kp2.octave];

            /*if(fabs(ratioDist-ratioOctave)>ratioFactor)
                continue;*/
            if(ratioDist*ratioFactor<ratioOctave || ratioDist>ratioOctave*ratioFactor)
                continue;

            // Triangulation is succesfull
            MapPoint* pMP = new MapPoint(x3D,mpCurrentKeyFrame,mpMap);

            pMP->AddObservation(mpCurrentKeyFrame,idx1);            
            pMP->AddObservation(pKF2,idx2);

            mpCurrentKeyFrame->AddMapPoint(pMP,idx1);
            pKF2->AddMapPoint(pMP,idx2);

            pMP->ComputeDistinctiveDescriptors();

            pMP->UpdateNormalAndDepth();

            mpMap->AddMapPoint(pMP);
            mlpRecentAddedMapPoints.push_back(pMP);

            nnew++;
        }
    }
}

void LocalMapping::SearchInNeighbors()
{
    // Retrieve neighbor keyframes
    int nn = 10;
    if(mbMonocular)
        nn=20;
    const vector<KeyFrame*> vpNeighKFs = mpCurrentKeyFrame->GetBestCovisibilityKeyFrames(nn);
    vector<KeyFrame*> vpTargetKFs;
    for(vector<KeyFrame*>::const_iterator vit=vpNeighKFs.begin(), vend=vpNeighKFs.end(); vit!=vend; vit++)
    {
        KeyFrame* pKFi = *vit;
        if(pKFi->isBad() || pKFi->mnFuseTargetForKF == mpCurrentKeyFrame->mnId)
            continue;
        vpTargetKFs.push_back(pKFi);
        pKFi->mnFuseTargetForKF = mpCurrentKeyFrame->mnId;

        // Extend to some second neighbors
        const vector<KeyFrame*> vpSecondNeighKFs = pKFi->GetBestCovisibilityKeyFrames(5);
        for(vector<KeyFrame*>::const_iterator vit2=vpSecondNeighKFs.begin(), vend2=vpSecondNeighKFs.end(); vit2!=vend2; vit2++)
        {
            KeyFrame* pKFi2 = *vit2;
            if(pKFi2->isBad() || pKFi2->mnFuseTargetForKF==mpCurrentKeyFrame->mnId || pKFi2->mnId==mpCurrentKeyFrame->mnId)
                continue;
            vpTargetKFs.push_back(pKFi2);
        }
    }


    // Search matches by projection from current KF in target KFs
    ORBmatcher matcher;
    vector<MapPoint*> vpMapPointMatches = mpCurrentKeyFrame->GetMapPointMatches();
    for(vector<KeyFrame*>::iterator vit=vpTargetKFs.begin(), vend=vpTargetKFs.end(); vit!=vend; vit++)
    {
        KeyFrame* pKFi = *vit;

        matcher.Fuse(pKFi,vpMapPointMatches);
    }

    // Search matches by projection from target KFs in current KF
    vector<MapPoint*> vpFuseCandidates;
    vpFuseCandidates.reserve(vpTargetKFs.size()*vpMapPointMatches.size());

    for(vector<KeyFrame*>::iterator vitKF=vpTargetKFs.begin(), vendKF=vpTargetKFs.end(); vitKF!=vendKF; vitKF++)
    {
        KeyFrame* pKFi = *vitKF;

        vector<MapPoint*> vpMapPointsKFi = pKFi->GetMapPointMatches();

        for(vector<MapPoint*>::iterator vitMP=vpMapPointsKFi.begin(), vendMP=vpMapPointsKFi.end(); vitMP!=vendMP; vitMP++)
        {
            MapPoint* pMP = *vitMP;
            if(!pMP)
                continue;
            if(pMP->isBad() || pMP->mnFuseCandidateForKF == mpCurrentKeyFrame->mnId)
                continue;
            pMP->mnFuseCandidateForKF = mpCurrentKeyFrame->mnId;
            vpFuseCandidates.push_back(pMP);
        }
    }

    matcher.Fuse(mpCurrentKeyFrame,vpFuseCandidates);


    // Update points
    vpMapPointMatches = mpCurrentKeyFrame->GetMapPointMatches();
    for(size_t i=0, iend=vpMapPointMatches.size(); i<iend; i++)
    {
        MapPoint* pMP=vpMapPointMatches[i];
        if(pMP)
        {
            if(!pMP->isBad())
            {
                pMP->ComputeDistinctiveDescriptors();
                pMP->UpdateNormalAndDepth();
            }
        }
    }

    // Update connections in covisibility graph
    mpCurrentKeyFrame->UpdateConnections();
}

cv::Mat LocalMapping::ComputeF12(KeyFrame *&pKF1, KeyFrame *&pKF2)
{
    cv::Mat R1w = pKF1->GetRotation();
    cv::Mat t1w = pKF1->GetTranslation();
    cv::Mat R2w = pKF2->GetRotation();
    cv::Mat t2w = pKF2->GetTranslation();

    cv::Mat R12 = R1w*R2w.t();
    cv::Mat t12 = -R1w*R2w.t()*t2w+t1w;

    cv::Mat t12x = SkewSymmetricMatrix(t12);

    const cv::Mat &K1 = pKF1->mK;
    const cv::Mat &K2 = pKF2->mK;


    return K1.t().inv()*t12x*R12*K2.inv();
}

void LocalMapping::RequestStop()
{
    unique_lock<mutex> lock(mMutexStop);
    mbStopRequested = true;
    unique_lock<mutex> lock2(mMutexNewKFs);
    mbAbortBA = true;
}

bool LocalMapping::Stop()
{
    unique_lock<mutex> lock(mMutexStop);
    if(mbStopRequested && !mbNotStop)
    {
        mbStopped = true;
        cout << "Local Mapping STOP" << endl;
        return true;
    }

    return false;
}

bool LocalMapping::isStopped()
{
    unique_lock<mutex> lock(mMutexStop);
    return mbStopped;
}

bool LocalMapping::stopRequested()
{
    unique_lock<mutex> lock(mMutexStop);
    return mbStopRequested;
}

void LocalMapping::Release()
{
    unique_lock<mutex> lock(mMutexStop);
    unique_lock<mutex> lock2(mMutexFinish);
    if(mbFinished)
        return;
    mbStopped = false;
    mbStopRequested = false;
    for(list<KeyFrame*>::iterator lit = mlNewKeyFrames.begin(), lend=mlNewKeyFrames.end(); lit!=lend; lit++)
        delete *lit;
    mlNewKeyFrames.clear();

    cout << "Local Mapping RELEASE" << endl;
}

bool LocalMapping::AcceptKeyFrames()
{
    unique_lock<mutex> lock(mMutexAccept);
    return mbAcceptKeyFrames;
}

void LocalMapping::SetAcceptKeyFrames(bool flag)
{
    unique_lock<mutex> lock(mMutexAccept);
    mbAcceptKeyFrames=flag;
}

bool LocalMapping::SetNotStop(bool flag)
{
    unique_lock<mutex> lock(mMutexStop);

    if(flag && mbStopped)
        return false;

    mbNotStop = flag;

    return true;
}

void LocalMapping::InterruptBA()
{
    mbAbortBA = true;
}

void LocalMapping::KeyFrameCulling()
{

    if(ConfigParam::GetRealTimeFlag())
    {
        if(GetFlagCopyInitKFs())
            return;
    }
    SetFlagCopyInitKFs(true);

    // Check redundant keyframes (only local keyframes)
    // A keyframe is considered redundant if the 90% of the MapPoints it sees, are seen
    // in at least other 3 keyframes (in the same or finer scale)
    // We only consider close stereo points
    vector<KeyFrame*> vpLocalKeyFrames = mpCurrentKeyFrame->GetVectorCovisibleKeyFrames();

    KeyFrame* pOldestLocalKF = mlLocalKeyFrames.front();
    KeyFrame* pPrevLocalKF = pOldestLocalKF->GetPrevKeyFrame();
    KeyFrame* pNewestLocalKF = mlLocalKeyFrames.back();
    // Test log
    if(pOldestLocalKF->isBad()) cerr<<"pOldestLocalKF is bad, check 1. id: "<<pOldestLocalKF->mnId<<endl;
    if(pPrevLocalKF) if(pPrevLocalKF->isBad()) cerr<<"pPrevLocalKF is bad, check 1. id: "<<pPrevLocalKF->mnId<<endl;
    if(pNewestLocalKF->isBad()) cerr<<"pNewestLocalKF is bad, check 1. id: "<<pNewestLocalKF->mnId<<endl;

    for(vector<KeyFrame*>::iterator vit=vpLocalKeyFrames.begin(), vend=vpLocalKeyFrames.end(); vit!=vend; vit++)
    {
        KeyFrame* pKF = *vit;
        if(pKF->mnId==0)
            continue;

        // Don't cull the oldest KF in LocalWindow,
        // And the KF before this KF
        if(pKF == pOldestLocalKF || pKF == pPrevLocalKF)
            continue;

        // Check time between Prev/Next Keyframe, if larger than 0.5s(for local)/3s(others), don't cull
        // Note, the KF just out of Local is similarly considered as Local
        KeyFrame* pPrevKF = pKF->GetPrevKeyFrame();
        KeyFrame* pNextKF = pKF->GetNextKeyFrame();
        if(pPrevKF && pNextKF && !GetVINSInited())
        {
            if(fabs(pNextKF->mTimeStamp - pPrevKF->mTimeStamp) > /*0.2*/0.5)
                continue;
        }
        // Don't drop the KF before current KF
        if(pKF->GetNextKeyFrame() == mpCurrentKeyFrame)
            continue;
        if(pKF->mTimeStamp >= mpCurrentKeyFrame->mTimeStamp - 0.11)
            continue;

        if(pPrevKF && pNextKF)
        {
            double timegap=0.51;
            if(GetVINSInited() && pKF->mTimeStamp < mpCurrentKeyFrame->mTimeStamp - 4.0)
                timegap = 3.01;

            if(fabs(pNextKF->mTimeStamp - pPrevKF->mTimeStamp) > timegap)
                continue;
        }

        const vector<MapPoint*> vpMapPoints = pKF->GetMapPointMatches();

        int nObs = 3;
        const int thObs=nObs;
        int nRedundantObservations=0;
        int nMPs=0;
        for(size_t i=0, iend=vpMapPoints.size(); i<iend; i++)
        {
            MapPoint* pMP = vpMapPoints[i];
            if(pMP)
            {
                if(!pMP->isBad())
                {
                    if(!mbMonocular)
                    {
                        if(pKF->mvDepth[i]>pKF->mThDepth || pKF->mvDepth[i]<0)
                            continue;
                    }

                    nMPs++;
                    if(pMP->Observations()>thObs)
                    {
                        const int &scaleLevel = pKF->mvKeysUn[i].octave;
                        const mapMapPointObs/*map<KeyFrame*, size_t>*/ observations = pMP->GetObservations();
                        int nObs=0;
                        for(mapMapPointObs/*map<KeyFrame*, size_t>*/::const_iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
                        {
                            KeyFrame* pKFi = mit->first;
                            if(pKFi==pKF)
                                continue;
                            const int &scaleLeveli = pKFi->mvKeysUn[mit->second].octave;

                            if(scaleLeveli<=scaleLevel+1)
                            {
                                nObs++;
                                if(nObs>=thObs)
                                    break;
                            }
                        }
                        if(nObs>=thObs)
                        {
                            nRedundantObservations++;
                        }
                    }
                }
            }
        }  

        if(nRedundantObservations>0.9*nMPs)
            pKF->SetBadFlag();
    }

    SetFlagCopyInitKFs(false);
}

cv::Mat LocalMapping::SkewSymmetricMatrix(const cv::Mat &v)
{
    return (cv::Mat_<float>(3,3) <<             0, -v.at<float>(2), v.at<float>(1),
            v.at<float>(2),               0,-v.at<float>(0),
            -v.at<float>(1),  v.at<float>(0),              0);
}

void LocalMapping::RequestReset()
{
    {
        unique_lock<mutex> lock(mMutexReset);
        mbResetRequested = true;
    }

    while(1)
    {
        {
            unique_lock<mutex> lock2(mMutexReset);
            if(!mbResetRequested)
                break;
        }
        usleep(3000);
    }
}

void LocalMapping::ResetIfRequested()
{
    unique_lock<mutex> lock(mMutexReset);
    if(mbResetRequested)
    {
        mlNewKeyFrames.clear();
        mlpRecentAddedMapPoints.clear();
        mbResetRequested=false;

        mlLocalKeyFrames.clear();

        // Add resetting init flags
        mbVINSInited = false;
        mbFirstTry = true;
    }
}

void LocalMapping::RequestFinish()
{
    unique_lock<mutex> lock(mMutexFinish);
    mbFinishRequested = true;
}

bool LocalMapping::CheckFinish()
{
    unique_lock<mutex> lock(mMutexFinish);
    return mbFinishRequested;
}

void LocalMapping::SetFinish()
{
    unique_lock<mutex> lock(mMutexFinish);
    mbFinished = true;    
    unique_lock<mutex> lock2(mMutexStop);
    mbStopped = true;
}

bool LocalMapping::isFinished()
{
    unique_lock<mutex> lock(mMutexFinish);
    return mbFinished;
}

} //namespace ORB_SLAM
