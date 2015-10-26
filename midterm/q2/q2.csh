#!/bin/csh

mv file.txt file.txt.old

foreach file ( *.txt )
   mv $file $file.old
end

foreach file ( *.txt )
   set root = `basename $file .txt`
   mv $root $root.dat

   set new = `echo $file | sed 's/txt/dat/'`
   mv $root $new
end

foreach file ( $* )
   mv $file $file.old
end
