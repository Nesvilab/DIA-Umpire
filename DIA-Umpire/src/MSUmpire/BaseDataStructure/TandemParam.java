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

import org.apache.commons.io.FilenameUtils;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class TandemParam extends DBSearchParam {
    
    public TandemParam(SearchInstrumentType type) {
        defaultType = type;
        SetParameter(type);
    }
    public void SetCombineFileName(String filename, String tag) {
        CombinedPepXML = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(filename) + "interact-" + FilenameUtils.getBaseName(filename) + tag + ".tandem.combine.pep.xml");
        CombinedProt = FilenameUtils.getFullPath(filename) + FilenameUtils.getBaseName(filename) + tag + ".tandem.Qcombine.prot.xml";
    }

    @Override
    public void SetResultFilePath(String mzXMLfile) {
        SpectrumPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(mzXMLfile) + FilenameUtils.getName(mzXMLfile));
        PepXMLPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(mzXMLfile) + FilenameUtils.getBaseName(mzXMLfile) + ".tandem.pep.xml");
        InteractPepXMLPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(mzXMLfile) + "interact-" + FilenameUtils.getBaseName(mzXMLfile) + ".tandem.pep.xml");
        ProtXMLPath = InteractPepXMLPath.replace(".pep.xml", ".prot.xml");
        parameterPath = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(mzXMLfile) + FilenameUtils.getBaseName(mzXMLfile) + ".tandem.param");
        RawSearchResult = FilenameUtils.separatorsToUnix(FilenameUtils.getFullPath(mzXMLfile) + FilenameUtils.getBaseName(mzXMLfile) + ".tandem");
    }
}
