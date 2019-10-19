# 3D RF Ultrasound Volume Reconstruction

The code provided was produced as part of the Ph.D thesis ["Statistical Image Processing of Medical Ultrasound Radio Frequency Data "](https://mediatum.ub.tum.de/1100919). The software has multiple purposes: First, to record B-mode or RF sequences in combination with tracking data. Second, to reconstruct 3D volumes from raw RF ultrasound data (radio-frequency data) or conventional B-mode. Third, there was some preliminary functionality to reconstruct flow velocity from Doppler ultrasound video streams. 
To acquire the data, the ultrasound transducer was tracked with an optical tracking system.

![3D Ultrasound Freehand System](https://github.com/TJKlein/3D_RFUltrasound_Reconstruction/blob/master/3DUSFreehand.png)

RF data is the more or less unprocessed signal at the beginning of processing pipeline that generates the commonly known ultrasound B-mode images. Given the comparably high resolution of the data and the raw signal characteristics make RF data attractive for statistical image processing.  


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
