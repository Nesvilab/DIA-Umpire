/*
 * This file is part of DIA-Umpire.
 *
 * DIA-Umpire is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the Free Software 
 * Foundation, either version 3 of the License, or (at your option) any later 
 * version.
 *
 * DIA-Umpire is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
 * more details.
 *
 * You should have received a copy of the GNU General Public License along 
 * with DIA-Umpire. If not, see <https://www.gnu.org/licenses/>.
 *
 */
package MSUmpire.DIA;

import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.LCMSPeakStructure.LCMSPeakDIAMS2;
import MSUmpire.LCMSPeakStructure.LCMSPeakMS1;
import MSUmpire.PSMDataStructure.LCMSID;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PeakDataStructure.PeakCluster;
import java.util.HashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.logging.log4j.LogManager;

/**
 * For a given isolation window, extract grouped fragments for all peptide ions in the isolation window 
 * @author Chih-Chiang Tsou
 */
public class DIA_window_Quant implements  Runnable{

    public LCMSPeakDIAMS2 DIAWindow;
    HashMap<Integer, Integer> ScanClusterMap_Q1;
    HashMap<Integer, Integer> ScanClusterMap_Q2;
    HashMap<Integer, String> ScanClusterMap_Q3;
    String Q1Name;
    String Q2Name;    
    String Q3Name;
    LCMSPeakMS1 ms1lcms;
    LCMSID IDsummary;
    int NoThread=1;
    
    public DIA_window_Quant(String Q1Name,String Q2Name,String Q3Name,HashMap<Integer, Integer> ScanClusterMap_Q1, HashMap<Integer, Integer> ScanClusterMap_Q2, HashMap<Integer, String> ScanClusterMap_Q3, LCMSPeakMS1 ms1lcms,LCMSPeakDIAMS2 DIAWindow,LCMSID IDsummary, int NoThreads){
        this.ScanClusterMap_Q1=ScanClusterMap_Q1;
        this.ScanClusterMap_Q2=ScanClusterMap_Q2;
        this.ScanClusterMap_Q3=ScanClusterMap_Q3;
        this.Q1Name=Q1Name;
        this.Q2Name=Q2Name;
        this.Q3Name=Q3Name;
        this.ms1lcms=ms1lcms;
        this.DIAWindow=DIAWindow;
        this.IDsummary=IDsummary;
        this.NoThread=NoThreads;
    }
    @Override
    public void run() {
       
        if(!DIAWindow.ReadPeakCluster()){
            LogManager.getRootLogger().error("Reading Peak cluster result for " + DIAWindow.ScanCollectionName + " failed");
            return;
        }
        ExecutorService executorPool;
        executorPool = Executors.newFixedThreadPool(NoThread);
        
        //For each identified peptide ion, extract the precursor feature and grouped fragments from the isolation window
        for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
            if (DIAWindow.DIA_MZ_Range.getX() <= pepIonID.GetPeakMz(2)&& DIAWindow.DIA_MZ_Range.getY() >= pepIonID.ObservedMz) {
                DIAMapClusterUnit mapunit = new DIAMapClusterUnit(pepIonID, Q1Name, Q2Name, Q3Name, ScanClusterMap_Q1, ScanClusterMap_Q2, ScanClusterMap_Q3, ms1lcms, DIAWindow);
                executorPool.execute(mapunit);
            }
        }
        executorPool.shutdown();

        try {
            executorPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            LogManager.getRootLogger().info("interrupted..");
        }
        if (DIAWindow.datattype != SpectralDataType.DataType.pSMART) {
            if (!DIAWindow.ReadPrecursorFragmentClu2Cur()) {
                LogManager.getRootLogger().error("Reading precursor-fragment results for " + DIAWindow.ScanCollectionName + " failed");
                return;
            }

            for (PepIonID pepIonID : IDsummary.GetPepIonList().values()) {
                for (PeakCluster cluster : pepIonID.MS1PeakClusters) {
                    if (DIAWindow.DIA_MZ_Range.getX() <= cluster.GetMaxMz() && DIAWindow.DIA_MZ_Range.getY() >= cluster.TargetMz()) {
                        DIAWindow.ExtractFragmentForPeakCluser(cluster);
                    }
                }
                for (PeakCluster ms2cluster : pepIonID.MS2UnfragPeakClusters) {
                    if (DIAWindow.DIA_MZ_Range.getX() <= ms2cluster.TargetMz() && DIAWindow.DIA_MZ_Range.getY() >= ms2cluster.TargetMz() && DIAWindow.PeakClusters.size()>=ms2cluster.Index) {
                        PeakCluster cluster = DIAWindow.PeakClusters.get(ms2cluster.Index - 1);
                        if (cluster.TargetMz() == ms2cluster.TargetMz() || cluster.Charge == ms2cluster.Charge) {
                            DIAWindow.ExtractFragmentForUnfragPeakCluser(cluster);
                        }
                    }
                }
            }
        }
        DIAWindow.ClearAllPeaks();
    }
    
}
