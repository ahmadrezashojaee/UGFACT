# UGFACT
Underground Gas Flow simulAtions with Coupled bio-geochemical reacTions (UGFACT)

This MATLAB-based framework couples reservoir simulation with bio-geochemical batch reactions using MRST and PHREEQC.

Developed by **Ahmadreza Shojaee** during his PhD at **Heriot-Watt University**, funded through the **James Watt Scholarship**.

## Acknowledgments

This work was conducted under the supervision of  
**Dr. Saeed Ghanbari**, **Dr. Gang Wang**, and **Prof. Eric Mackay**  
at the **Institute of GeoEnergy Engineering**, School of Energy, Geoscience, Infrastructure and Society (EGIS), Heriot-Watt University.

Their guidance and support throughout the development of this framework are gratefully acknowledged.

---

## 💻 System Requirements

- **MATLAB** (R2022a or later recommended)
- **MRST 2024b** from SINTEF
- **PHREEQC / IPhreeqc** from USGS
- Windows OS (required for IPhreeqc COM integration)

---

## 📦 Installation Instructions

### 1. Download MRST
Download MRST from SINTEF:
https://www.sintef.no/projectweb/mrst/download/

### 2. Framework Compatibility
This framework is tested with **MRST-2024b**. It may work with other versions, but this has not been verified. Please test carefully if using other versions.

### 3. Integrate UGFACT Module
Download or clone the `UGFACT` module and place it into the MRST modules folder:
```
...\mrst-2024b\modules\UGFACT
```

### 4. Download IPhreeqc (COM version)
From the USGS website:
https://water.usgs.gov/water-resources/software/PHREEQC/index.html

Recommended version:
https://water.usgs.gov/water-resources/software/PHREEQC/IPhreeqcCOM-3.8.6-17100-x64.msi

### 5. Install IPhreeqc
Run the `.msi` installer and install IPhreeqc on your system.

### 6. Replace IPhreeqc DLLs in MRST
Navigate to the installed IPhreeqc location:
```
C:\Program Files\USGS\IPhreeqcCOM 3.8.6-17100\bin\
```

Copy the following files:
- `IPhreeqcCOM.dll`
- `IPhreeqcCOM.tlb`

Paste and replace the versions located in:
```
...\mrst-2024b\modules\UGFACT\Modified_Models\
```

### 7. Run MRST startup
Run the MRST startup script from MATLAB:
```matlab
run('...\mrst-2024b\startup.m')
```

### 8. Enable Parallel Execution (Optional)
To enable parallel batch calculations using PHREEQC:
```matlab
model.parpool = true;
```

### 9. Run Examples
Navigate to the examples directory:
```
...\mrst-2024b\modules\UGFACT\examples
```

Run any example script to verify the installation.

---

## 🧾 Licensing & Credits

- **PHREEQC**: Developed by the USGS – [USGS Software License](https://www.usgs.gov/products/software/overview)
- **MRST**: Developed by SINTEF – [GPLv3 License](https://www.gnu.org/licenses/gpl-3.0.html)

Please cite PHREEQC and MRST appropriately in academic or published work.

---

© 2025 Ahmadreza Shojaee – Heriot-Watt University
