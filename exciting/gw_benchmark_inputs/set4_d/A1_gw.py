# TODO Add mixbasis, freqgrid, barecoul, selfenergy
gw_string_template = \
"""  
   <gw
    taskname="{taskname}"
    nempty="{nempty}"
    ngridq="{ngridq}"
    skipgnd="{skipgnd}"
    >
  
    <mixbasis
      lmaxmb="4"
      epsmb="1.d-3"
      gmb="1.d0"
    ></mixbasis>
  
    <freqgrid
      nomeg="32"
      freqmax="1.0"
    ></freqgrid>
  
    <barecoul
      pwm="2.0"
      stctol="1.d-16"
      barcevtol="0.1"
    ></barecoul>
  
    <selfenergy
      actype="pade"
      singularity="mpb"
    ></selfenergy>
  
   </gw>
"""