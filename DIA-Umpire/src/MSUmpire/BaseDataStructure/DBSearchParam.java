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
package MSUmpire.BaseDataStructure;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public abstract class DBSearchParam implements Cloneable{

    public float FragPPM;
    public float PrecursorPPM;
    public int MinNoPeaksScoring;
    public int MinNoPeaks;
    public int TotalPeaks;
    public String SpectrumPath;
    public String RawSearchResult;
    public String InteractPepXMLPath;
    public String PepXMLPath;
    public String ProtXMLPath;
    public String CombinedPepXML;
    public String CombinedProt;
    public String FastaPath;
    public String DecoyFasta;
    public String OutputSeqPath;
    public int MissCleavage = 1;
    public boolean SemiCleavage = false;
    public boolean NonSpecificCleavage=false;
    public boolean IsotopeError = false;    
    public String parameterPath;
    public String templateParamFile;
    public SearchInstrumentType defaultType;
    public int NoCPUs = 2;
    public float PepFDR = 0.01f;
    public float ProtFDR = 0.01f;
    public boolean Overwrite = false;
    public String DecoyPrefix="rev_";
    
    public enum SearchInstrumentType {
        Orbitrap,
        TOF5600,
        QExactive,
        };

    @Override
    public Object clone() throws CloneNotSupportedException {
        return super.clone();
    }

    public abstract void SetResultFilePath(String mzXMLfile);
    public abstract void SetCombineFileName(String filename, String tag);
            
    protected void SetParameter(SearchInstrumentType type) {
        switch (type) {
            case QExactive: {
                PrecursorPPM = 10;
                FragPPM = 20;
                MinNoPeaksScoring = 3;
                MinNoPeaks = 15;
                TotalPeaks = 140;
                MissCleavage = 1;
                break;
            }
            case Orbitrap: {
                PrecursorPPM = 10;
                FragPPM = 500;
                MinNoPeaksScoring = 3;
                MinNoPeaks = 15;
                TotalPeaks = 100;
                break;
            }

            case TOF5600: {
                PrecursorPPM = 30;
                FragPPM = 40;
                MinNoPeaksScoring = 3;
                MinNoPeaks = 3;
                TotalPeaks = 140;
                MissCleavage = 1;
                SemiCleavage = false;
                break;
            }
        }
    }
}
