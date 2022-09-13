function [] = dc_make_movie(indir,files,frameRate,outname)
    
    olddir = pwd;
    cd(indir)
    cmd = sprintf('E:\\Software\\Misc\\ffmpeg\\bin\\ffmpeg.exe -b:v 25000000 -q:v 1 -r %d -f image2 -i %s -q:v 1 -g 1 %s',frameRate,files,outname);
    system(cmd);
    cd(olddir);