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
import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeakDataStructure.PrecursorFragmentPairEdge;
import MSUmpire.PeakDataStructure.SortedCorrFrag;
import com.compomics.util.experiment.biology.ions.ElementaryIon;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

/**
 * Preprocessing to generate pseudo MS/MS spectra
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PseudoMSMSProcessing implements Runnable {

    InstrumentParameter parameter;
    ArrayList<PrecursorFragmentPairEdge> fragments;
    public PeakCluster Precursorcluster;
    float growth = 1f;

    public PseudoMSMSProcessing(PeakCluster ms1cluster, InstrumentParameter parameter) {
        this.parameter = parameter;
        this.fragments = ms1cluster.GroupedFragmentPeaks;
        this.Precursorcluster = ms1cluster;
    }

    public void DeisotopingForPeakClusterFragment() {
        ArrayList<PrecursorFragmentPairEdge> newfragments = new ArrayList<>();
        boolean[] fragmentmarked = new boolean[fragments.size()];
        Arrays.fill(fragmentmarked, Boolean.TRUE);
        PrecursorFragmentPairEdge currentmaxfragment = fragments.get(0);
        int currentmaxindex = 0;
        for (int i = 1; i < fragments.size(); i++) {
            if (InstrumentParameter.CalcPPM(fragments.get(i).FragmentMz, currentmaxfragment.FragmentMz) > parameter.MS2PPM) {
                fragmentmarked[currentmaxindex] = false;
                currentmaxindex = i;
                currentmaxfragment = fragments.get(i);
            } else if (fragments.get(i).Intensity > currentmaxfragment.Intensity) {
                currentmaxindex = i;
                currentmaxfragment = fragments.get(i);
            }
        }
        fragmentmarked[currentmaxindex] = false;
        for (int i = 0; i < fragments.size(); i++) {
            if (!fragmentmarked[i]) {
                fragmentmarked[i] = true;
                PrecursorFragmentPairEdge startfrag = fragments.get(i);

                boolean groupped = false;
                for (int charge = 2; charge >= 1; charge--) {
                    float lastint = startfrag.Intensity;
                    boolean found = false;
                    for (int pkidx = 1; pkidx < 5; pkidx++) {
                        float targetmz = startfrag.FragmentMz + (float) pkidx / charge;
                        for (int j = i + 1; j < fragments.size(); j++) {
                            if (!fragmentmarked[j]) {
                                PrecursorFragmentPairEdge targetfrag = fragments.get(j);
                                if (InstrumentParameter.CalcPPM(targetfrag.FragmentMz, targetmz) < parameter.MS2PPM * (pkidx * 0.5 + 1)) {
                                    if (targetfrag.Intensity < lastint) {
                                        fragmentmarked[j] = true;
                                        lastint = targetfrag.Intensity;
                                        found = true;
                                        break;
                                    }
                                } else if (targetfrag.FragmentMz > targetmz) {
                                    break;
                                }
                            }
                        }
                        if (!found) {
                            break;
                        }
                    }
                    if (found) {
                        groupped = true;
                        //convert to charge 1 m/z
                        startfrag.FragmentMz = startfrag.FragmentMz * charge - (charge - 1) * (float) ElementaryIon.proton.getTheoreticMass();
                        if (startfrag.FragmentMz <= Precursorcluster.NeutralMass()) {
                            newfragments.add(startfrag);
                        }
                    }
                }
                if (!groupped) {
                    newfragments.add(startfrag);
                }
            }
        }

        fragments = newfragments;
        SortFragmentByMZ();
    }

    public void RemoveMatchedFrag(HashMap<Integer, ArrayList<PrecursorFragmentPairEdge>> MatchedFragmentMap) {
        ArrayList<PrecursorFragmentPairEdge> newlist = new ArrayList<>();
        for (PrecursorFragmentPairEdge fragmentClusterUnit : fragments) {
            if (!MatchedFragmentMap.containsKey(fragmentClusterUnit.PeakCurveIndexB)) {
                newlist.add(fragmentClusterUnit);
            }
        }
        fragments = newlist;
    }

    public void BoostComplementaryIon() {
        float totalmass = (float) (Precursorcluster.TargetMz() * Precursorcluster.Charge - Precursorcluster.Charge * ElementaryIon.proton.getTheoreticMass());
        boolean[] fragmentmarked = new boolean[fragments.size()];
        Arrays.fill(fragmentmarked, Boolean.FALSE);
        for (int i = 0; i < fragments.size(); i++) {
            PrecursorFragmentPairEdge fragmentClusterUnit = fragments.get(i);
            if (!fragmentmarked[i]) {
                fragmentmarked[i] = true;
                ArrayList<PrecursorFragmentPairEdge> GroupedFragments = new ArrayList<>();
                GroupedFragments.add(fragmentClusterUnit);
                float complefrag1 = (float) (totalmass - fragmentClusterUnit.FragmentMz + 2f * ElementaryIon.proton.getTheoreticMass());
                if (complefrag1 >= fragmentClusterUnit.FragmentMz) {
                    for (int j = i + 1; j < fragments.size(); j++) {
                        if (!fragmentmarked[j]) {
                            PrecursorFragmentPairEdge fragmentClusterUnit2 = fragments.get(j);
                            if (InstrumentParameter.CalcPPM(complefrag1, fragmentClusterUnit2.FragmentMz) < parameter.MS2PPM) {
                                GroupedFragments.add(fragmentClusterUnit2);
                                fragmentmarked[j] = true;
                            } else if (fragmentClusterUnit2.FragmentMz > complefrag1) {
                                break;
                            }
                        }
                    }
                }

                if (GroupedFragments.size() > 1) {
                    PrecursorFragmentPairEdge bestfragment = GroupedFragments.get(0);
                    for (PrecursorFragmentPairEdge fragment : GroupedFragments) {
                        if (fragment.Intensity > bestfragment.Intensity) {
                            bestfragment = fragment;
                        }
                    }
                    for (PrecursorFragmentPairEdge fragment : GroupedFragments) {
                        fragment.ComplementaryFragment = true;
                        fragment.Intensity = bestfragment.Intensity * growth;
                        fragment.Correlation = bestfragment.Correlation;
                        fragment.ApexDelta = bestfragment.ApexDelta;
                        fragment.RTOverlapP = bestfragment.RTOverlapP;
                        fragment.FragmentMS1Rank = bestfragment.FragmentMS1Rank;
                        fragment.FragmentMS1RankScore = bestfragment.FragmentMS1RankScore;
                    }
                }
            }
        }
    }

    public void IdentifyComplementaryIon(float totalmass) {
        boolean[] fragmentmarked = new boolean[fragments.size()];
        Arrays.fill(fragmentmarked, Boolean.FALSE);
        for (int i = 0; i < fragments.size(); i++) {
            PrecursorFragmentPairEdge fragmentClusterUnit = fragments.get(i);
            if (!fragmentmarked[i]) {
                fragmentmarked[i] = true;
                ArrayList<PrecursorFragmentPairEdge> GroupedFragments = new ArrayList<>();
                GroupedFragments.add(fragmentClusterUnit);
                float complefrag1 = (float) (totalmass - fragmentClusterUnit.FragmentMz + 2f * ElementaryIon.proton.getTheoreticMass());
                if (complefrag1 >= fragmentClusterUnit.FragmentMz) {
                    for (int j = i + 1; j < fragments.size(); j++) {
                        if (!fragmentmarked[j]) {
                            PrecursorFragmentPairEdge fragmentClusterUnit2 = fragments.get(j);
                            if (InstrumentParameter.CalcPPM(complefrag1, fragmentClusterUnit2.FragmentMz) < parameter.MS2PPM) {
                                GroupedFragments.add(fragmentClusterUnit2);
                                fragmentmarked[j] = true;
                            } else if (fragmentClusterUnit2.FragmentMz > complefrag1) {
                                break;
                            }
                        }
                    }
                }

                for (PrecursorFragmentPairEdge fragment : GroupedFragments) {
                    fragment.ComplementaryFragment = true;
                }
            }
        }
    }

    public XYPointCollection GetScan() {
        XYPointCollection Scan = new XYPointCollection();
        for (PrecursorFragmentPairEdge fragmentClusterUnit : fragments) {
            if (parameter.AdjustFragIntensity) {
                Scan.AddPointKeepMaxIfCloseValueExisted(fragmentClusterUnit.FragmentMz, fragmentClusterUnit.Intensity * fragmentClusterUnit.Correlation * fragmentClusterUnit.Correlation, parameter.MS2PPM);                
            }
            else{
                Scan.AddPointKeepMaxIfCloseValueExisted(fragmentClusterUnit.FragmentMz, fragmentClusterUnit.Intensity, parameter.MS2PPM);
            }
        }
        return Scan;
    }

    private void SortFragmentByMZ() {
        Collections.sort(fragments, new Comparator<PrecursorFragmentPairEdge>() {
            @Override
            public int compare(PrecursorFragmentPairEdge o1, PrecursorFragmentPairEdge o2) {
                return Float.compare(o1.FragmentMz, o2.FragmentMz);
            }
        });
    }

    @Override
    public void run() {
        if (fragments.size() > 1) {
            SortFragmentByMZ();
            if (parameter.BoostComplementaryIon) {
                DeisotopingForPeakClusterFragment();
                BoostComplementaryIon();
            }
        }
    }
}
