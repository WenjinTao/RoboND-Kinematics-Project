#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *
import numpy as np


def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
        def rot_x(angle):
            return Matrix([[ 1,         0,            0],
                           [ 0, cos(angle), -sin(angle)],
                           [ 0, sin(angle),  cos(angle)]])
        def rot_y(angle):
            return Matrix([[ cos(angle),  0, sin(angle)],
                           [          0,  1,          0],
                           [-sin(angle),  0, cos(angle)]])
        def rot_z(angle):
            return Matrix([[  cos(angle), -sin(angle), 0],
                           [  sin(angle),  cos(angle), 0],
                           [      0,        0,         1]])

        def DH_matrix(alpha, a, d, q):
            return Matrix([[            cos(q),           -sin(q),           0,             a],
                           [ sin(q)*cos(alpha), cos(q)*cos(alpha), -sin(alpha), -sin(alpha)*d],
                           [ sin(q)*sin(alpha), cos(q)*sin(alpha),  cos(alpha),  cos(alpha)*d],
                           [                 0,                 0,           0,             1]])

        # Initialize service response
        joint_trajectory_list = []

        # Define DH param symbols
        alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7') # twist angle
        a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7') # joint length
        d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8') # joint offset
                    
        # Joint angle symbols
        q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8') # joint variables: theta_i
        roll, pitch, yaw = symbols('roll pitch yaw')
        
#        theta1,theta2,theta3,theta4,theta5,theta6 = 0,0,0,0,0,0

        # Modified DH params
        s = {alpha0:     0, a0:      0, d1:  0.75,
             alpha1: -pi/2, a1:   0.35, d2:     0, q2: q2-pi/2,
             alpha2:     0, a2:   1.25, d3:     0,
             alpha3: -pi/2, a3: -0.054, d4:  1.50,
             alpha4:  pi/2, a4:      0, d5:     0,
             alpha5: -pi/2, a5:      0, d6:     0,
             alpha6:     0, a6:      0, d7: 0.303, q7: 0}

        # Define Modified DH Transformation matrix
        
        T0_1 = DH_matrix(alpha0, a0, d1, q1).subs(s)        
        T1_2 = DH_matrix(alpha1, a1, d2, q2).subs(s)        
        T2_3 = DH_matrix(alpha2, a2, d3, q3).subs(s)
        
#        T3_4 = DH_matrix(alpha3, a3, d4, q4).subs(s)        
#        T4_5 = DH_matrix(alpha4, a4, d5, q5).subs(s)        
#        T5_6 = DH_matrix(alpha5, a5, d6, q6).subs(s)        
#        T6_G = DH_matrix(alpha6, a6, d7, q7).subs(s)
        
        # Rrpy = Homogeneous RPY rotation between base_link and gripper_link
        # The overall RPY (Roll Pitch Yaw) rotation between base_link and
        # gripper_link must be equal to the product of individual rotations
        # between respective links: R0_6 = Rrpy

        R_roll = rot_x(roll)
        R_pitch = rot_y(pitch)
        R_yaw = rot_z(yaw)
       
        # Correction to account of orientation difference between definition of gripper link in URDF versus DH convention
        # [yaw-pitch-roll to rotation matirx](http://planning.cs.uiuc.edu/node102.html)
        R_z = rot_z(pi)
        R_y = rot_y(-pi/2)
        R_corr = R_z * R_y

        Rrpy = R_yaw * R_pitch * R_roll * R_corr

        T0_3 = T0_1 * T1_2 * T2_3
#        T0_G = T0_3 * T3_4 * T4_5 * T5_6 * T6_G
#        T_total = T0_G * R_corr

        R0_3 = T0_3[0:3, 0:3]
        
        R3_6 = R0_3**(-1) * Rrpy
        
#        print 'got-R3_6'

        for i in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

            # Extract end-effector position and orientation from request
            # px,py,pz = end-effector position
            # roll, pitch, yaw = end-effector orientation
            px = req.poses[i].position.x
            py = req.poses[i].position.y
            pz = req.poses[i].position.z

            (ee_roll, ee_pitch, ee_yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[i].orientation.x, req.poses[i].orientation.y,
                    req.poses[i].orientation.z, req.poses[i].orientation.w])
     
            # Calculate joint angles using Geometric IK method
            
            # 1. calculate WC: (WC_x, WC_y, WC_z)
            P_EE = Matrix([px, py, pz])
            P_WC = P_EE - Rrpy.evalf(subs={roll: ee_roll, pitch: ee_pitch, yaw: ee_yaw}) * Matrix([0, 0, s[d7]])
            
            P_WC_x, P_WC_y, P_WC_z = P_WC[0], P_WC[1], P_WC[2]

            # 2. calculate theta1:
            theta1 = atan2(P_WC_y, P_WC_x).evalf()
            
            # 3. calculate theta2:
            # S1: distance between the x-y-plane projections of 2 and 5
            # s2: P_WC_Z - d1
            # s3: distance between 2, 5
            # s4: distance between 3, 5
            # beta1: atan2(s2, s1)
            # beta2: angle_325
            # beta3: angle_235
            # beta4: pi/2 - angle_23(5^o)
            
            S1 = sqrt(P_WC_x**2 + P_WC_y**2) - s[a1]
            S2 = P_WC_z - s[d1]
            S3 = sqrt(S1**2 + S2**2)
            S4 = sqrt(s[d4]**2 + s[a3]**2)
            beta1 = atan2(S2, S1)
            cos_beta2 = (s[a2]**2 + S3**2 - S4**2)/(2*s[a2]*S3)
#            print cos_beta2
            sin_beta2 = sqrt(1-cos_beta2**2)
            beta2 = atan2(sin_beta2, cos_beta2)
            
            theta2 = (pi/2 - beta1 - beta2).evalf()             

            # 4. calculate theta3:
            cos_beta3 = (s[a2]**2 + S4**2 - S3**2)/(2*s[a2]*S4)
            sin_beta3 = sqrt(1 - cos_beta3**2)
            beta3 = atan2(sin_beta3, cos_beta3)
            beta4 = atan2(-s[a3], s[d4]) # s[a3] is negative
            
            theta3 = (pi/2 - beta4 - beta3).evalf()

            # 5. calculate theta4, 5, 6:
            R3_6_num = R3_6.evalf(subs={q1: theta1, q2: theta2, q3: theta3, roll: ee_roll, pitch: ee_pitch, yaw: ee_yaw})
            
            r13 = R3_6_num[0,2]
            r33 = R3_6_num[2,2]
            
            r21 = R3_6_num[1,0]
            r22 = R3_6_num[1,1]
            
            r23 = R3_6_num[1,2]
            
            theta4 = atan2(r33, -r13).evalf()
            theta6 = atan2(-r22, r21).evalf()
            theta5 = atan2(sqrt(1 - r23**2), r23).evalf()
                        

            print ('px, py, pz: ', px, py, pz)
            print ('roll, pitch, yaw: ', ee_roll, ee_pitch, ee_yaw)            
            print ('p_wc_x, p_wc_y, p_wc_z: ', P_WC_x, P_WC_y, P_WC_z)
            print ('theta1~6: ', theta1, theta2, theta3, theta4, theta5, theta6)


            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
            joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
            joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
