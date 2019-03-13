#!/bin/bash
source ../RunAtStart.sh
source ../SourceYoshidaSan.sh
source ../RunAtStart.sh
./bin/Linux-g++/WCSim WCSimHKHybridmPMT_14374_50Hz_FullTank$1MeV.mac tuningNominal.mac > /disk01/usr5/bquilain/tt$1.txt
