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

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.LCMSPeakStructure.LCMSPeakMS1;
import MSUmpire.PSMDataStructure.FragmentPeak;
import MSUmpire.PSMDataStructure.PepIonID;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeakDataStructure.PrecursorFragmentPairEdge;
import com.compomics.util.experiment.biology.Ion;
import com.compomics.util.experiment.biology.ions.ElementaryIon;
import com.compomics.util.experiment.biology.ions.PeptideFragmentIon;
import org.apache.commons.io.FilenameUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Thread unit for determining  peak cluster matched fragments for identified peptide ion for quantification
 * 
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class DIAAssignQuantUnit implements Runnable {

    PepIonID pepIonID;
    LCMSPeakMS1 ms1lcms;
    double protonMass = ElementaryIon.proton.getTheoreticMass();
    InstrumentParameter parameter;

    public DIAAssignQuantUnit(PepIonID pepIonID, LCMSPeakMS1 ms1lcms, InstrumentParameter parameter) {
        this.pepIonID = pepIonID;
        this.ms1lcms = ms1lcms;
        this.parameter = parameter;
    }

    @Override
    public void run() {        
        PeakCluster targetCluster = null;
        
        //Get highest intensity peak cluster
        for (PeakCluster peakCluster : pepIonID.MS1PeakClusters) {
            if (targetCluster == null || targetCluster.PeakHeight[0] < peakCluster.PeakHeight[0]) {
                targetCluster = peakCluster;
            }
        }        
        pepIonID.CreateQuantInstance(ms1lcms.MaxNoPeakCluster);

        if (targetCluster != null) {
            pepIonID.PeakArea = targetCluster.PeakArea;
            pepIonID.PeakHeight = targetCluster.PeakHeight;
            pepIonID.PeakClusterScore = targetCluster.MS1Score;
            pepIonID.PeakRT = targetCluster.PeakHeightRT[0];
            pepIonID.ObservedMz = targetCluster.mz[0];
        } else {
            //if no MS1 peak cluster found, use MS2 unfragmented peak cluster instead
            for (PeakCluster peakCluster : pepIonID.MS2UnfragPeakClusters) {
                if (targetCluster == null || targetCluster.PeakHeight[0] < peakCluster.PeakHeight[0]) {
                    targetCluster = peakCluster;
                }
            }
            if (targetCluster != null) {
                pepIonID.PeakRT = targetCluster.PeakHeightRT[0];
                pepIonID.ObservedMz = targetCluster.mz[0];
            }
        }
        //pepIonID.PeakClusterIndex = targetCluster.Index;

        if (targetCluster!=null && ms1lcms.datattype != SpectralDataType.DataType.pSMART) {
            MatchFragmentByTargetCluster(targetCluster);
            pepIonID.RemoveRedundantFrag();
            if (pepIonID.FragmentPeaks.isEmpty() && Math.max(pepIonID.MaxProbability, pepIonID.TargetedProbability()) > 0.8f) {
                LogManager.getRootLogger().warn("Warning: " + pepIonID.ModSequence + "(MaxProb: " + pepIonID.MaxProbability + ") does not have matched fragment in "+FilenameUtils.getBaseName(ms1lcms.ParentmzXMLName));
                //MatchFragment();
            }
            pepIonID.ClearPepFragFactory();
        }
    }

    //Determine matched fragments
    private void MatchFragmentByTargetCluster(PeakCluster peakCluster) {
        for (Ion frag : pepIonID.GetFragments()) {
            PrecursorFragmentPairEdge bestfragment = null;
            //Singly charged framgnet ion
            float targetmz = (float) frag.getTheoreticMz(1);
            for (PrecursorFragmentPairEdge fragmentClusterUnit : peakCluster.GroupedFragmentPeaks) {
                if (InstrumentParameter.CalcPPM(targetmz, fragmentClusterUnit.FragmentMz) <= parameter.MS2PPM) {
                    if (bestfragment == null || fragmentClusterUnit.Correlation > bestfragment.Correlation) {
                        bestfragment = fragmentClusterUnit;
                    }
                }
            }
           
            //Singly-charged fragment ion is found
            if (bestfragment != null) {
                FragmentPeak fragmentpeak = new FragmentPeak();
                fragmentpeak.ObservedMZ = bestfragment.FragmentMz;
                fragmentpeak.FragMZ = targetmz;
                fragmentpeak.corr = bestfragment.Correlation;
                fragmentpeak.intensity = bestfragment.Intensity;
                fragmentpeak.ApexDelta = bestfragment.ApexDelta;
                fragmentpeak.RTOverlapP = bestfragment.RTOverlapP;
                fragmentpeak.Charge = 1;
                fragmentpeak.ppm = InstrumentParameter.CalcSignedPPM(bestfragment.FragmentMz, targetmz);
                fragmentpeak.IonType = frag.getSubTypeAsString() + ((PeptideFragmentIon) frag).getNumber();
                pepIonID.FragmentPeaks.add(fragmentpeak);
            }

            //doubly charged fragment ion
            targetmz = (float) frag.getTheoreticMz(2);
            bestfragment = null;
            for (PrecursorFragmentPairEdge fragmentClusterUnit : peakCluster.GroupedFragmentPeaks) {
                if (InstrumentParameter.CalcPPM(targetmz, fragmentClusterUnit.FragmentMz) <= parameter.MS2PPM) {
                    if (bestfragment == null || fragmentClusterUnit.Correlation > bestfragment.Correlation) {
                        bestfragment = fragmentClusterUnit;
                    }
                }
            }

            //Doubly-charged fragment ion is found
            if (bestfragment != null) {
                FragmentPeak fragmentpeak = new FragmentPeak();
                fragmentpeak.ObservedMZ = bestfragment.FragmentMz;
                fragmentpeak.FragMZ = targetmz;
                fragmentpeak.corr = bestfragment.Correlation;
                fragmentpeak.intensity = bestfragment.Intensity;
                fragmentpeak.ApexDelta = bestfragment.ApexDelta;
                fragmentpeak.RTOverlapP = bestfragment.RTOverlapP;
                fragmentpeak.Charge = 2;
                fragmentpeak.ppm = InstrumentParameter.CalcSignedPPM(bestfragment.FragmentMz, targetmz);
                fragmentpeak.IonType = frag.getSubTypeAsString() + ((PeptideFragmentIon) frag).getNumber();
                pepIonID.FragmentPeaks.add(fragmentpeak);
            }
        }
    }    
}
