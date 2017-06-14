# ZEUS-MP v2.0

This repository contains the web page and source code for ZEUS-MP 2.  This code was originally written by Dr. John Hayes at the Laboratory for Computational Astrophysics at the University of Illinois at Urbana-Champaign (until 2000) and the University of California at San Diego (2000 onward), and is based on the numerical methods developed for previous versions of the ZEUS code.  The code capabilities and numerical methods are described in [Hayes et al. (1996), ApJ Supplement, 165, 188-228](http://adsabs.harvard.edu/abs/2006ApJS..165..188H).

This repository is being maintained by Brian O'Shea ([oshea@msu.edu](mailto:oshea@msu.edu)) as a service to the community, since the server hosting the original ZEUS-MP 2 website is defunct and the code is not being actively developed.  Questions about the codebase should be directed to John Hayes ([jchayes63@gmail.com](jchayes63@gmail.com)), the primary author of ZEUS-MP 2.  

This repository contains two primary directories:

```CODE/``` contains the source code for ZEUS-MP 2.  This includes a subdirectory ```zeusmp2```, which is the last developed version of the code.  A second subdirectory, ```released_vers```, contains gzipped tar files corresponding to numbered versions of the code, which contain bug fixes and some new features (most prominently HDF5 support in Version 2.1.0).  Version 2.0.0 was released on November 16, 2005, and the following versions were released over roughly the next 16 months, culiminating in Version 2.1.2, which was released on January 25, 2007.  Each tarball contains a filed named RELEASE_NOTES, which explains the bug fixes and new features.

```DOCS/``` contains the entire website for ZEUS-MP 2, which includes documentation for Version 2.1.2 (which is accurate for previous versions as well), including html and PDF versions of the documentation.  The file [DOCS/MAN/intro.html](DOCS/MAN/intro.html) may be particularly useful to the new user of ZEUS-MP 2, as it explains the relationship between this codebase and other versions of the ZEUS code (starting with ZEUS-2D, [Stone and Norman (1992)](http://adsabs.harvard.edu/abs/1992ApJS...80..753S))

### LICENSE

ZEUS-MP Version 2 is released to the community under the University of Illinois/NCSA Open Source License (code written prior to 2000) and the University of California/BSD License (code written after 2000). These licenses are compatible with each other.  Please consult the documentation for full versions of the license if you have any questions.

### WARNING

ZEUS-MP 2 is distributed **as-is**, with no guarantees regarding user support or bug fixes.
