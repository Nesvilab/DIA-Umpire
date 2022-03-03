/* 
 * Author: Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 *             Nesvizhskii Lab, Department of Computational Medicine and Bioinformatics, 
 *             University of Michigan, Ann Arbor
 *
 * Copyright 2014 University of Michigan, Ann Arbor, MI
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package MSUmpire.DIA;

import MSUmpire.LCMSPeakStructure.LCMSPeakDIAMS2;
import MSUmpire.LCMSPeakStructure.LCMSPeakMS1;
import MSUmpire.PSMDataStructure.PSM;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PeakDataStructure.PeakCluster;
import java.util.ArrayList;
import java.util.HashMap;
import org.apache.commons.io.FilenameUtils;
import org.apache.logging.log4j.LogManager;

/**
 * Thread unit for assigning MS1 peak cluster and matched MS2 fragment peak for identified peptide ion
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class DIAMapClusterUnit implements Runnable{

    PepIonID pepIonID;
    HashMap<Integer, Integer> ScanClusterMap_Q1;
    HashMap<Integer, Integer> ScanClusterMap_Q2;
    HashMap<Integer, String> ScanClusterMap_Q3;
    String Q1Name;
    String Q2Name;    
    String Q3Name;
    LCMSPeakMS1 ms1lcms;
    LCMSPeakDIAMS2 DIAWindow;
    
    public  DIAMapClusterUnit(PepIonID pepIonID, String Q1Name,String Q2Name,String Q3Name,HashMap<Integer, Integer> ScanClusterMap_Q1, HashMap<Integer, Integer> ScanClusterMap_Q2, HashMap<Integer, String> ScanClusterMap_Q3, LCMSPeakMS1 ms1lcms, LCMSPeakDIAMS2 DIAWindow){
        this.pepIonID=pepIonID;
        this.ScanClusterMap_Q1=ScanClusterMap_Q1;
        this.ScanClusterMap_Q2=ScanClusterMap_Q2;
        this.ScanClusterMap_Q3=ScanClusterMap_Q3;
        this.Q1Name=Q1Name;
        this.Q2Name=Q2Name;
        this.Q3Name=Q3Name;
        this.ms1lcms=ms1lcms;
        this.DIAWindow=DIAWindow;
    }
    @Override
    public void run() {
        //For each identified PSM
         for (PSM psm : pepIonID.GetPSMList()) {
            int ClusterIndex = -1;
            if (psm.GetRawNameString() == null ? Q1Name == null : psm.GetRawNameString().equals(FilenameUtils.getBaseName(Q1Name))) {
                if(!ScanClusterMap_Q1.containsKey(psm.ScanNo)){
                   LogManager.getRootLogger().error("ScanClusterMapping error");
                   LogManager.getRootLogger().error("ScanClusterMapping "+Q1Name+" doesn't have "+ psm.SpecNumber);
                    System.exit(3);
                }
                //Get cluster index fro Q1
                ClusterIndex = ScanClusterMap_Q1.get(psm.ScanNo);
                PeakCluster Cluster = ms1lcms.PeakClusters.get(ClusterIndex - 1);
                Cluster.Identified=true;
                if (!pepIonID.MS1PeakClusters.contains(Cluster)) {
                    pepIonID.MS1PeakClusters.add(Cluster);
                }
            } else if (psm.GetRawNameString() == null ? Q2Name == null : psm.GetRawNameString().equals(FilenameUtils.getBaseName(Q2Name))) {
                if(!ScanClusterMap_Q2.containsKey(psm.ScanNo)){
                    LogManager.getRootLogger().error("ScanClusterMapping error");
                    LogManager.getRootLogger().error("ScanClusterMapping "+Q2Name+" doesn't have "+ psm.SpecNumber);
                    System.exit(3);
                }
                
                //Get cluster index fro Q2
                ClusterIndex = ScanClusterMap_Q2.get(psm.ScanNo);
                PeakCluster Cluster = ms1lcms.PeakClusters.get(ClusterIndex - 1);
                Cluster.Identified=true;
                if (!pepIonID.MS1PeakClusters.contains(Cluster)) {
                    pepIonID.MS1PeakClusters.add(Cluster);
                }
             } else if (psm.GetRawNameString() == null ? Q3Name == null : psm.GetRawNameString().equals(FilenameUtils.getBaseName(Q3Name))) {
                 String WindowClusterIndex = ScanClusterMap_Q3.get(psm.ScanNo);
                 if(!ScanClusterMap_Q3.containsKey(psm.ScanNo)){
                    LogManager.getRootLogger().error("ScanClusterMapping error");
                    LogManager.getRootLogger().error("ScanClusterMapping "+Q3Name+" doesn't have "+ psm.SpecNumber);
                    System.exit(3);
                }
                 if (WindowClusterIndex.split(";").length == 2) {
                     String windowname = WindowClusterIndex.split(";")[0];
                     if (windowname.equals(DIAWindow.WindowID)) {
                         //Get cluster index fro Q3
                         ClusterIndex = Integer.parseInt(WindowClusterIndex.split(";")[1]);
                         PeakCluster Cluster = DIAWindow.PeakClusters.get(ClusterIndex - 1);
                         Cluster.Identified = true;
                         pepIonID.MS2UnfragPeakClusters.add(Cluster);
                         ArrayList<PeakCluster> ms1list = ms1lcms.FindPeakClustersByMassRTRange(Cluster.NeutralMass(), Cluster.Charge, Cluster.startRT, Cluster.endRT);
                         for (PeakCluster ms1cluster : ms1list) {
                             ms1cluster.Identified = true;
                             if (!pepIonID.MS1PeakClusters.contains(ms1cluster)) {
                                 pepIonID.MS1PeakClusters.add(ms1cluster);
                             }
                         }
                     }
                 } else {
                     //Get cluster index fro Q3
                     ClusterIndex = Integer.parseInt(WindowClusterIndex);
                     if (DIAWindow.UnFragIonClu2Cur.containsKey(ClusterIndex)) {
                         PeakCluster Cluster = DIAWindow.PeakClusters.get(ClusterIndex - 1);
                         if (Cluster.Charge == psm.Charge && Math.abs(Cluster.TargetMz() - psm.ObserPrecursorMz()) < 0.01f && Math.abs(Cluster.PeakHeightRT[0] - psm.RetentionTime) < 0.1f) {
                             Cluster.Identified = true;
                             pepIonID.MS2UnfragPeakClusters.add(Cluster);
                         }
                         ArrayList<PeakCluster> ms1list = ms1lcms.FindPeakClustersByMassRTRange(Cluster.NeutralMass(), Cluster.Charge, Cluster.startRT, Cluster.endRT);
                         for (PeakCluster ms1cluster : ms1list) {
                             ms1cluster.Identified = true;
                             if (!pepIonID.MS1PeakClusters.contains(ms1cluster)) {
                                 pepIonID.MS1PeakClusters.add(ms1cluster);
                             }
                         }
                     }
                 }                 
             }
        }       
         
         if(pepIonID.MS1PeakClusters.isEmpty() && pepIonID.MS2UnfragPeakClusters.isEmpty()){
             LogManager.getRootLogger().trace("Cannot find feature for identified peptide ion : "+pepIonID.GetKey()+" mz: "+pepIonID.ObservedMz+" in isolation window "+DIAWindow.WindowID);
         }
    }

}
