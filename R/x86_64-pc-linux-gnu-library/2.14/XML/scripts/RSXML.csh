if(`test -n "-L/usr/lib/x86_64-linux-gnu -lxml2"`) then

if(${?LD_LIBRARY_PATH}) then
    setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:-L/usr/lib/x86_64-linux-gnu -lxml2
else
   setenv LD_LIBRARY_PATH -L/usr/lib/x86_64-linux-gnu -lxml2
endif

endif