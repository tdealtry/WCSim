#!/bin/bash
#./bin/Linux-g++/WCSim WCSimHKHybridmPMT14374100HzFullTank3MeV.mac tuningNominal.mac > /disk01/usr5/bquilain/tt3.txt
#./bin/Linux-g++/WCSim WCSimHKHybridmPMT14374100HzFullTank4MeV.mac tuningNominal.mac > /disk01/usr5/bquilain/tt4.txt
#./bin/Linux-g++/WCSim WCSimHKHybridmPMT14374100HzFullTank5MeV.mac tuningNominal.mac > /disk01/usr5/bquilain/tt5.txt
#./bin/Linux-g++/WCSim WCSimHKHybridmPMT14374100HzFullTank6MeV.mac tuningNominal.mac > /disk01/usr5/bquilain/tt6.txt
#./bin/Linux-g++/WCSim WCSimHKHybridmPMT14374100HzFullTank8MeV.mac tuningNominal.mac > /disk01/usr5/bquilain/tt8.txt
#./bin/Linux-g++/WCSim WCSimHKHybridmPMT14374100HzFullTank10MeV.mac tuningNominal.mac > /disk01/usr5/bquilain/tt10.txt
#./bin/Linux-g++/WCSim WCSimHKHybridmPMT14374100HzFullTank15MeV.mac tuningNominal.mac > /disk01/usr5/bquilain/tt15.txt

#./bin/Linux-g++/WCSim WCSimHKHybridmPMT14374100HzFullTank3MeV_noDN.mac tuningNominal.mac > /disk01/usr5/bquilain/tt3_nodn.txt
#./bin/Linux-g++/WCSim WCSimHKHybridmPMT14374100HzFullTank4MeV_noDN.mac tuningNominal.mac > /disk01/usr5/bquilain/tt4_nodn.txt
#./bin/Linux-g++/WCSim WCSimHKHybridmPMT14374100HzFullTank5MeV_noDN.mac tuningNominal.mac > /disk01/usr5/bquilain/tt5_nodn.txt
#./bin/Linux-g++/WCSim WCSimHKHybridmPMT14374100HzFullTank6MeV_noDN.mac tuningNominal.mac > /disk01/usr5/bquilain/tt6_nodn.txt
#./bin/Linux-g++/WCSim WCSimHKHybridmPMT14374100HzFullTank8MeV_noDN.mac tuningNominal.mac > /disk01/usr5/bquilain/tt8_nodn.txt
#./bin/Linux-g++/WCSim WCSimHKHybridmPMT14374100HzFullTank10MeV_noDN.mac tuningNominal.mac > /disk01/usr5/bquilain/tt10_nodn.txt
#./bin/Linux-g++/WCSim WCSimHKHybridmPMT14374100HzFullTank15MeV_noDN.mac tuningNominal.mac > /disk01/usr5/bquilain/tt15_nodn.txt
for i in 3 4 5 6 8 10 15
#for i in 3
do
#    rm -f WCSimHKHybridmPMT10PC14374100HzFullTank${i}MeV.mac
    cp WCSimHKHybridmPMT14374100HzFullTank${i}MeV.mac temp.mac
#    sed 's/HyperK_HybridmPMT/HyperK_HybridmPMT10PC/g' temp.mac > temp2.mac
#    sed 's/HyperK_HybridmPMT/HyperK/g' temp.mac > temp2.mac
    sed 's/hkhybridmpmt14374100Hz/hkhybridmpmt1437450Hz/g' temp.mac > WCSimHKHybridmPMT_14374_50Hz_FullTank${i}MeV.mac
#    cp WCSimHKHybridmPMT14374100HzFullTank${i}MeV.mac WCSimHKHybridmPMT10PC14374100HzFullTank${i}MeV.mac
#    sed 's/HyperK_HybridmPMT/HyperK_HybridmPMT10PC/g' WCSimHKHybridmPMT10PC14374100HzFullTank${i}MeV.mac
#    sed 's/hkhybridmpmt/hkhybridmpmt10pc/g' WCSimHKHybridmPMT10PC14374100HzFullTank${i}MeV.mac
    
done	 
