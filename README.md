# 3D RF Ultrasound Volume Reconstruction

The code provided in this repository was largely produced as part of the Ph.D thesis ["Statistical Image Processing of Medical Ultrasound Radio Frequency Data "](https://mediatum.ub.tum.de/1100919). The software has multiple purposes: First, to record B-mode or RF sequences in combination with tracking data. Second, to reconstruct 3D volumes from raw RF ultrasound data (radio-frequency data) or conventional B-mode. Third, there was some preliminary functionality to reconstruct flow velocity from Doppler ultrasound video streams. 
To acquire the data, the ultrasound transducer was tracked with an optical tracking system.

![3D Ultrasound Freehand System](https://github.com/TJKlein/3D_RFUltrasound_Reconstruction/blob/master/3DUSFreehand.png)

  
## Background

Conventional ultrasound images, commonly referred to as B-Mode, are the result of many processing
steps optimizing data for visual assessment by physicians. However, at the core of ultrasound imaging pipeline lies the radio frequency (RF) data. Just lately, RF data has become more readily available to the research community such that its potential has not fully unveiled yet. From a data processing standpoint using RF data over B-Mode suggests many advantages. First of all, it is generally much richer in information due to the comparably higher resolution. Furthermore, it is not affected by non-linear post-processing steps such as log-compression and proprietary filter algorithms that change the noise statistics for reasons of improved visual appeal. In addition, it has nice probabilistic properties facilitating various ways of distributional modeling of ultrasound specific texture patterns, referred to as speckle noise.

![RF to Bmode pipeline](https://github.com/TJKlein/3D_RFUltrasound_Reconstruction/blob/master/RFtoBmode.png)




## Citation
If you use this code or find it somehow useful for your research, I would appreciate citation:


```
@inproceedings{Klein2012StatisticalIP,
  title={Statistical Image Processing of Medical Ultrasound Radio Frequency Data},
  author={Tassilo Klein},
  year={2012}
}
```

```
@InProceedings{10.1007/978-3-642-33415-3_52,
author="Klein, T.
and Hansson, M.
and Navab, Nassir",
editor="Ayache, Nicholas
and Delingette, Herv{\'e}
and Golland, Polina
and Mori, Kensaku",
title="Modeling of Multi-View 3D Freehand Radio Frequency Ultrasound",
booktitle="Medical Image Computing and Computer-Assisted Intervention -- MICCAI 2012",
year="2012",
publisher="Springer Berlin Heidelberg",
address="Berlin, Heidelberg",
pages="422--429"
}
```

```
@InProceedings{10.1007/978-3-319-24571-3_71,
author="Klein, Tassilo
and Wells, William M.",
editor="Navab, Nassir
and Hornegger, Joachim
and Wells, William M.
and Frangi, Alejandro",
title="RF Ultrasound Distribution-Based Confidence Maps",
booktitle="Medical Image Computing and Computer-Assisted Intervention -- MICCAI 2015",
year="2015",
publisher="Springer International Publishing",
address="Cham",
pages="595--602"
}
```

```
@InProceedings{10.1007/978-3-642-33415-3_52,
author="Klein, T.
and Hansson, M.
and Navab, Nassir",
editor="Ayache, Nicholas
and Delingette, Herv{\'e}
and Golland, Polina
and Mori, Kensaku",
title="Modeling of Multi-View 3D Freehand Radio Frequency Ultrasound",
booktitle="Medical Image Computing and Computer-Assisted Intervention -- MICCAI 2012",
year="2012",
publisher="Springer Berlin Heidelberg",
address="Berlin, Heidelberg",
pages="422--429"
}
```
