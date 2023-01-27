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

package MSUmpire.SpectrumParser;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.ScanData;
import uk.ac.ebi.jmzml.model.mzml.BinaryDataArray;
import uk.ac.ebi.jmzml.model.mzml.CVParam;
import uk.ac.ebi.jmzml.model.mzml.ParamGroup;
import uk.ac.ebi.jmzml.model.mzml.Precursor;
import uk.ac.ebi.jmzml.model.mzml.PrecursorList;
import uk.ac.ebi.jmzml.model.mzml.ScanList;
import uk.ac.ebi.jmzml.model.mzml.Spectrum;

/**
 * Thread unit to covert mzML spectrum to ScanData class
 * @author Chih-Chiang Tsou
 */
public class mzMLSpecConverter implements Runnable{

    public ScanData spec=null;
    Spectrum jmzMLSpec=null;
    InstrumentParameter parameter=null;
    public mzMLSpecConverter(Spectrum jmzMLSpec, InstrumentParameter parameter){
        this.jmzMLSpec=jmzMLSpec;
        this.parameter=parameter;
    }
    @Override
    public void run() {
        getScanFromJMzMLSpec();
        spec.Preprocessing(parameter);
    }
    
    public void getScanFromJMzMLSpec() {
        spec = new ScanData();

        String id = jmzMLSpec.getId();
        spec.MGFTitle = id;
        // SpecIndex
        spec.ScanNum=jmzMLSpec.getIndex();

        // scan number
        String[] idToken = id.split("\\s+");
        if (idToken.length > 0 && idToken[idToken.length - 1].matches("scan=\\d+")) {
            int scanNum = Integer.parseInt(idToken[idToken.length - 1].substring(5));
            spec.ScanNum = scanNum;
        }

        // MS Level
        CVParam msLevelParam = null;
        Boolean isCentroided = false;
        for (CVParam cvParam : jmzMLSpec.getCvParam()) {
            if (cvParam.getAccession().equals("MS:1000511")) // MS level
            {
                msLevelParam = cvParam;
            } else if (cvParam.getAccession().equals("MS:1000127")) // centroid spectrum
            {
                isCentroided = true;
            } else if (cvParam.getAccession().equals("MS:1000128")) // profile spectrum
            {
                isCentroided = false;
            }
        }
        spec.centroided = isCentroided;

        float RT=-1f;
        
        // Scan list to get monoisotopic m/z
         ScanList scanList = jmzMLSpec.getScanList();
        if (scanList != null && scanList.getScan().size() > 0) {            
            for (CVParam param : scanList.getScan().get(0).getCvParam()) {
                if (param.getAccession().equals("MS:1000016")) // retention time
                {
                    RT = Float.parseFloat(param.getValue());
                    if (param.getUnitName().equals("second")) {
                        RT = RT / 60f;
                    }
                }
            }
        }

        int msLevel = msLevelParam != null ? Integer.parseInt(msLevelParam.getValue()) : 0;
        spec.MsLevel = msLevel;
        spec.RetentionTime=RT;

        // Precursor
        float precursorMz = -1f;
        PrecursorList precursorList = jmzMLSpec.getPrecursorList();
        if (precursorList != null && precursorList.getCount().intValue() > 0 && precursorList.getPrecursor().get(0).getSelectedIonList() != null) {
            Precursor precursor = precursorList.getPrecursor().get(0);	// consider only the first precursor

            ParamGroup isolationWindowParams = precursor.getIsolationWindow();
            if (isolationWindowParams != null && isolationWindowParams.getCvParam() != null) {
                float isolationWindowTargetMz = 0f;
                float Loffset=0f;
                float Roffset=0f;
                for (CVParam param : isolationWindowParams.getCvParam()) {
                    if (param.getAccession().equals("MS:1000827")) // selected ion m/z
                    {
                        isolationWindowTargetMz = Float.parseFloat(param.getValue());	// assume that unit is m/z (MS:1000040)
                    }
                    if (param.getAccession().equals("MS:1000828")) //lower offset
                    {
                        Loffset = Float.parseFloat(param.getValue());	
                    }
                     if (param.getAccession().equals("MS:1000829")) // upper offset
                    {
                        Roffset = Float.parseFloat(param.getValue());	
                    }
                     spec.isolationWindowTargetMz = isolationWindowTargetMz;
                     spec.isolationWindowLoffset=Loffset;
                     spec.isolationWindowRoffset=Roffset;
                }
            }

            // precursor mz, charge
            int precursorCharge = 0;
            float precursorIntensity = 0;

            ParamGroup paramGroup = precursor.getSelectedIonList().getSelectedIon().get(0);
            for (CVParam param : paramGroup.getCvParam()) {
                if (precursorMz < 0 && param.getAccession().equals("MS:1000744")) // selected ion m/z
                {
                    precursorMz = Float.parseFloat(param.getValue());	// assume that unit is m/z (MS:1000040)
                } else if (param.getAccession().equals("MS:1000041")) // charge state
                {
                    precursorCharge = Integer.parseInt(param.getValue());
                } else if (param.getAccession().equals("MS:1000042")) // peak intensity
                {
                    precursorIntensity = Float.parseFloat(param.getValue());
                }	//MS:1000511
            }
            spec.PrecursorMz = precursorMz;
            spec.PrecursorIntensity = precursorIntensity;
            spec.PrecursorCharge = precursorCharge;
            if (spec.isolationWindowTargetMz == 0f) {
                spec.isolationWindowTargetMz = precursorMz;
            }

            // activation method
            ParamGroup actMethodParams = precursor.getActivation();
            for (CVParam param : actMethodParams.getCvParam()) {
                spec.ActivationMethod = param.getValue();
                //System.out.println(param.getAccession());
                //System.out.println(param.getValue());
            }
        }

        // Peak list
        BinaryDataArray mzArray = null, intenArray = null;

        if (jmzMLSpec.getBinaryDataArrayList() != null && jmzMLSpec.getBinaryDataArrayList().getBinaryDataArray() != null) {
            for (BinaryDataArray array : jmzMLSpec.getBinaryDataArrayList().getBinaryDataArray()) {
                if (array.getEncodedLength() == 0) {
                    continue;
                }
                // check the cvParams
                for (CVParam param : array.getCvParam()) {
                    if (param.getAccession().equals("MS:1000514")) {
                        mzArray = array;
                        break;
                    }
                    if (param.getAccession().equals("MS:1000515")) {
                        intenArray = array;
                        break;
                    }
                }
                if (mzArray != null && intenArray != null) {
                    break;
                }
            }
        }

        if (mzArray != null && intenArray != null) {
            Number mzNumbers[] = mzArray.getBinaryDataAsNumberArray();
            Number intenNumbers[] = intenArray.getBinaryDataAsNumberArray();

            if (mzNumbers.length != intenNumbers.length) {
                System.err.println("Different sizes for m/z and intensity value arrays for spectrum" + jmzMLSpec.getId());
                System.exit(-1);
            }

            for (int i = 0; i < mzNumbers.length; i++) {
                spec.AddPoint(mzNumbers[i].floatValue(), intenNumbers[i].floatValue());
            }
        }             
    }
}
