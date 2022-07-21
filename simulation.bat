/* Comment: This is a batch file that will call and run the Mplus program, mplus_code.inp. Save it as “simulation.bat” (without the quotes). In order to run this file in Windows, under Start select Run and navigate to this file. You will need to edit the lines referring to C:\research\chapter to point to the directory in which you have saved the mplus_code.inp file. Prior to running the batch file, delete these comment lines. */

C:
cd\program files\mplus
mplus.exe C:\research\chapter\mplus_code.inp
move C:\Users\epspeg\AppData\Local\VirtualStore\Progra~1\Mplus\mplus_code.out C:\research\chapter
exit
