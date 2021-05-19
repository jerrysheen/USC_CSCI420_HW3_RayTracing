Assignment #3: Ray tracing

FULL NAME: TAO SHEN

MANDATORY FEATURES
------------------

<Under "Status" please indicate whether it has been implemented and is
functioning correctly.  If not, please explain the current status.>

Feature:                                 Status: finish? (yes/no)
-------------------------------------    -------------------------
1) Ray tracing triangles                  yes

2) Ray tracing sphere                     yes

3) Triangle Phong Shading               yes

4) Sphere Phong Shading                 yes

5) Shadows rays                          yes

6) Still images                          yes
   
7) Extra Credit (up to 20 points)
   !!! explain your extra credit here, if applicable !!!

Usage:  <input scenefile> [enable antialising y] [enable reflction y] [enable softshadow y] [output jpegname]\n
try enablt softshadow but not yet.

Antialising:  collect the pixel around the center ray(x + 0.5, y + 0.5) to (0.25, 0.25), (0.25, 0.75), (0.75, 0.75)(0.75, 0.25)
Reflection: using localphong_color * (1- fact) + ks * reflect_color * fact, to reduce the reflection effect..
Animation: ball motion.
