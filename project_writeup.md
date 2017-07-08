# Robotic arm: Pick & Place

## 1. Introduction

This project is to implement forward and inverse kinematics to control a Kuka KR210 robotic arm following a planned path. As shown in the following figure, the robot needs to pick a cylinder from the shelf and then place it into the bin. The moving trajectory is provided as a series of way points. To follow those points, we need to implement forward and inverse kinematics to obtain each joint angle. Then those angles are used to control this robotic arm.  

![][image3]

---

[//]: # "Image References"

[image1]: ./misc_images/frame_urdf.png
[image2]: ./misc_images/frame_dh.png
[image3]: ./misc_images/misc3.png
[figure_theta123]: ./misc_images/figure_theta123.png
[r3_6]: ./misc_images/r3_6.png


## 2. Forward Kinematic Analysis
### 2.1 DH Parameters Derivation 

The robot information is described as `URDF` files. The following figure shows the reference frame assignments and link dimensions in `URDF` file, from where we can derive the DH parameters.

![][image1]

Following the DH conventions and lecture notes, the DH frames are defined as follows.

![][image2]

By substituting the link dimensions from `URDF` file, the DH parameters can be derived:

| $i$  | $\alpha_{i-1} $ | $a_{i-1}$ | $d_i$ |  $q_i$  |
| :--: | :-------------: | :-------: | :---: | :-----: |
|  1   |        0        |     0     | 0.75  |         |
|  2   |      -pi/2      |   0.35    |   0   | q2-pi/2 |
|  3   |       *0*       |   1.25    |   0   |         |
|  4   |      -pi/2      |  -0.054   | 1.50  |         |
|  5   |      pi/2       |     0     |   0   |         |
|  6   |      -pi/2      |     0     |   0   |         |
|  G   |        0        |     0     | 0.303 |    0    |

Those DH parameters are defined as a dictionary:

```python
s = {alpha0:     0, a0:      0, d1:  0.75,
     alpha1: -pi/2, a1:   0.35, d2:     0, q2: q2-pi/2,
     alpha2:     0, a2:   1.25, d3:     0,
     alpha3: -pi/2, a3: -0.054, d4:  1.50,
     alpha4:  pi/2, a4:      0, d5:     0,
     alpha5: -pi/2, a5:      0, d6:     0,
     alpha6:     0, a6:      0, d7: 0.303, q7: 0}
```



#### 2. Transformation Matrices

A `DH_matrix(alpha, a, d, q)` function is defined as follows to create individual transformation matrices about each joint.

```python
def DH_matrix(alpha, a, d, q):
    return Matrix([[            cos(q),           -sin(q),           0,             a],
                   [ sin(q)*cos(alpha), cos(q)*cos(alpha), -sin(alpha), -sin(alpha)*d],
                   [ sin(q)*sin(alpha), cos(q)*sin(alpha),  cos(alpha),  cos(alpha)*d],
                   [                 0,                 0,           0,             1]])
```

Then the DH transformation matrices are implemented:

```python 
T0_1 = DH_matrix(alpha0, a0, d1, q1).subs(s)        
T1_2 = DH_matrix(alpha1, a1, d2, q2).subs(s)        
T2_3 = DH_matrix(alpha2, a2, d3, q3).subs(s)        
T3_4 = DH_matrix(alpha3, a3, d4, q4).subs(s)        
T4_5 = DH_matrix(alpha4, a4, d5, q5).subs(s)        
T5_6 = DH_matrix(alpha5, a5, d6, q6).subs(s)        
T6_G = DH_matrix(alpha6, a6, d7, q7).subs(s)

# Forward kinematics
T0_3 = T0_1 * T1_2 * T2_3
T0_G = T0_3 * T3_4 * T4_5 * T5_6 * T6_G
T_total = T0_G * R_corr
```

As the end-effector pose is provided, a homogeneous RPY(roll-pitch-yaw) rotation matrix can be constructed. Note that we need the follow the yaw-pitch-roll sequence when creating the matrix. `R_corr` is to account of orientation difference between the definition of gripper link in URDF versus DH convention.

```python
Rrpy = R_yaw * R_pitch * R_roll * R_corr
```



## 3. Inverse Kinematics Analysis 

### 3.1 Inverse Position Kinematics

First the location of the spherical wrist center (WC) can be calculated:

$$ \begin{bmatrix} P\_WC_x\\ P\_WC_y\\  P\_WC_z\end{bmatrix} = \begin{bmatrix} P\_EE_x\\ P\_EE_y\\  P\_EE_z\end{bmatrix}-Rrpy\begin{bmatrix} 0\\ 0\\ d_G\end{bmatrix}$$

where $P\_EE$ is the end-effector location. 

![alt text][figure_theta123]

Then the first 3 joint angles can be solved using trigonometry. As shown in the above figure:

$\theta_1 = arctan(P\_WC_y/P\_WC_x)$

$\theta_2=\pi/2-\beta_1-\beta_2$

$\theta_3=\pi/2-\beta_3-\beta_4$

where,

```python
# calculate theta1
theta1 = atan2(P_WC_y, P_WC_x).evalf()

# calculate theta2
S1 = sqrt(P_WC_x**2 + P_WC_y**2) - s[a1]
S2 = P_WC_z - s[d1]
S3 = sqrt(S1**2 + S2**2)
S4 = sqrt(s[d4]**2 + s[a3]**2)
beta1 = atan2(S2, S1)
cos_beta2 = (s[a2]**2 + S3**2 - S4**2)/(2*s[a2]*S3)
sin_beta2 = sqrt(1-cos_beta2**2)
beta2 = atan2(sin_beta2, cos_beta2)

theta2 = (pi/2 - beta1 - beta2).evalf()             

# calculate theta3:
cos_beta3 = (s[a2]**2 + S4**2 - S3**2)/(2*s[a2]*S4)
sin_beta3 = sqrt(1 - cos_beta3**2)
beta3 = atan2(sin_beta3, cos_beta3)
beta4 = atan2(-s[a3], s[d4]) # s[a3] is negative

theta3 = (pi/2 - beta4 - beta3).evalf()
```

### 3.2 Inverse Orientation Kinematics

After we have the first 3 joint angles, the rotation matrix `R3_6` of the last three joints can be calculated numerically:

```python
R0_3 = T0_3[0:3, 0:3]
R3_6 = R0_3**(-1) * Rrpy

R3_6_num = R3_6.evalf(subs={q1: theta1, q2: theta2, q3: theta3, roll: ee_roll, pitch: ee_pitch, yaw: ee_yaw})
```

While we can get its expression by:

`R3_6_exp = simplify(T3_4 * T4_5 * T5_6)[0:3, 0:3]`

which returns

![][r3_6]

Then the last three joint angles are

```python
r13 = R3_6_num[0,2]
r33 = R3_6_num[2,2]

r21 = R3_6_num[1,0]
r22 = R3_6_num[1,1]

r23 = R3_6_num[1,2]

theta4 = atan2(r33, -r13).evalf()
theta6 = atan2(-r22, r21).evalf()
theta5 = atan2(sqrt(1 - r23**2), r23).evalf()
```



## 4. Project Implementation

All the above methods are implemented in the `IK_server.py` file which is properly commented. The code is able to guide the robot to successfully complete the pick and place cycles. There are still several aspects need to be improved:

- Low efficiency of the **Calculating Inverse Kinematics** step. This may come from the symbolic computing part.
- The robotic arm is not moving fluently to its target. I'm not sure whether it is from the joint angle calculation or the simulator.
- The last 3 joints have some unnecessary rotations near the starting points under some situations. I think this is from IK calculation of  $\theta_{1,2,3}$ and should be solved by adding some angle limitations.
- Got this error {`Received 0 eef-poses from the plan No valid poses received`} once and failed that cycle. I think it came from the eef-poses sender.


