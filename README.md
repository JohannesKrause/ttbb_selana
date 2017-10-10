# ttbb_selana
The file *SELANA.C*  implements an experimental filter which removes overlap between tt+jets / ttbb from tt+jets events. 
The file has to be placed in the same directory as the run card.

compile instructions

1. set variable *SHERPA_PREFIX* to your installation path, e.g. 
    
            SHERPA_PREFIX=/home/name/software/Sherpa/rel-2-2-2
            
2. compile with 
    
            g++ -shared -g -I`$SHERPA_PREFIX/bin/Sherpa-config --incdir`  `$SHERPA_PREFIX/bin/Sherpa-config --ldflags`  -fPIC -o libSELANA.so SELANA.C
            

In order to use this filter, the created libary has to be linked to Sherpa and activated.
This can be done by including

  SHERPA_LDADD=SELANA;
  ANALYSIS SELANA; 

in the (run){  }(run) section of your run card.
If you want to use additional analysis modules, *SELANA* has to be the first one. 
In combination with Rivet this would mean   

ANALYSIS SELANA,Rivet;
           
