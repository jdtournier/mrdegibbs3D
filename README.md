This implements the 3D extension to the unringing method proposed by Kellner et al. (2016) [1]. This was presented by Thea Bautista at ISMRM 2021 [2]. 

_Note that this is a preview version. Its functionality is in the process of being merged into [the `mrdegibbs` command](https://mrtrix.readthedocs.io/en/latest/reference/commands/mrdegibbs.html) that currently implements the 2D version in [the main MRtrix3 repo](https://github.com/MRtrix3/mrtrix3)._

### Installation

You first need to have a source installation of [the main MRtrix3 repo](https://github.com/MRtrix3/mrtrix3) available. If you don't already, these are the minimum steps required:
```
git clone https://github.com/mrtrix3/mrtrix3.git
cd mrtrix3
./configure
```
You can add the `-nogui` option to the `./configure` call to avoid the checks for OpenGL & Qt if you don't intend to compile `mrview` (or `shview`).

You can then clone & compile this repo:
```
git clone https://github.com/jdtournier/mrdegibbs3D
cd mrdegibbs3D/
../build
```
For the `build` command above to work, it needs to invoke the main MRtrix3 `build` script (wherever this may be located) with the `mrdegibbs3D` repo as the current working directory. See [the MRtrix3 documentation on modules for full details](https://mrtrix.readthedocs.io/en/latest/tips_and_tricks/external_modules.html). 

If everything works as it should, you'll find the `degibbs3D` executable in the `bin/` folder. You can invoke it directly using its full path (e.g. `~/mrtrix3/mrdegibbs3D/bin/degibbs3D`), or add it to your `PATH`. 

**Note:** with the instructions above, you _cannot_ relocate the executable elsewhere since it will need to link to the main MRtrix3 library, and expects to find it in the same relative location as when it was compiled. You can however relocate it if you also ensure the main `libmrtrix3.so` library can still be found at the same _relative_ path to the executable itself. Alternatively, you can add the `-noshared` flag to the `./configure` call above to avoid the generation of the MRtrix3 library, which will result in a larger executable, but with fewer restrictions on relocation. 

---

1. [E. Kellner, B. Dhital, V. G. Kiselev, and M. Reisert, ‘Gibbs-ringing artifact removal based on local subvoxel-shifts’, Magnetic Resonance in Medicine, vol. 76, no. 5, pp. 1574–1581, 2016](https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.26054).
2. [T. Bautista, J. O’Muircheartaigh, J. V. Hajnal, and J.-D. Tournier, ‘Removal of Gibbs ringing artefacts for 3D acquisitions using subvoxel shifts’, in Proc. Intl. Soc. Mag. Reson. Med., online, May 2021, vol. 29, p. 3535](https://index.mirasmart.com/ISMRM2021/PDFfiles/3535.html).
