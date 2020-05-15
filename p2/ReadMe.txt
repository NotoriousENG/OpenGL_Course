**************
*Dependencies*
**************
1) Depends on tinyobjloader to load .obj files exported with Blender 2.79 and 2.8. 
Install with vcpkg, modify the path for CMake if needed.

***********
*DEMO GIFS*
***********
1) Camera: https://raw.githubusercontent.com/NotoriousENG/ogl_tutorial/master/KeyBoardInteraction.gif?token=AHWX43DEOX7BKSGRPHMRNAC6QQDQC
2) Keyboard Interaction: https://raw.githubusercontent.com/NotoriousENG/ogl_tutorial/master/Jump.gif?token=AHWX43HCJWDZPMLLQRDVIG26QQDQI
3) Teleporting: https://raw.githubusercontent.com/NotoriousENG/ogl_tutorial/master/Jump.gif?token=AHWX43HCJWDZPMLLQRDVIG26QQDQI

*********
*Issues *
*********
1) Pen Interaction: pressing shift and spinning pen works fine, other rotations do not work properly, not sure why. 
Probably has something to do with how I made the obj file.

2) Jumping: Since my pen's pivot is wrong so is the jumping, ball shoots from the origin.
