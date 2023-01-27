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
import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.BaseDataStructure.SpectralDataType;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.TreeMap;

/**
 * MALDI parser, parse a collection of pkl files (Data from Phil Andrews's lab)
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class MALDIDataParser {

    public ScanCollection scanCollection = null;
    public String filename;
    public InstrumentParameter parameter;
    public int TotalScan;
    public SpectralDataType.DataType datatype;
    public float cycletime=0.2f;
    public MALDIDataParser(String filename) throws IOException {
        parameter = new InstrumentParameter(InstrumentParameter.InstrumentType.Q_TOF);
        this.filename = filename;
        scanCollection = new ScanCollection(parameter.Resolution);
        scanCollection.Filename = filename;
        this.datatype = SpectralDataType.DataType.MALDI;                
    }

    public void Parse() throws FileNotFoundException, IOException {
        TreeMap<String, ScanData> SortedMap=new TreeMap<>();
        File folder =new File(filename);
                
        File[] listOfFiles = folder.listFiles();
        if (listOfFiles == null) {
            return;
        }
        int ScanNo=1;
        
        for (int i = 0; i < listOfFiles.length; i++) {
            if (listOfFiles[i].isFile()) {
                String files = listOfFiles[i].getName();
                if (files.toLowerCase().endsWith(".mgf")) {
                    MGFParser mgf = new MGFParser(folder + "/" + files);
                    for (Integer scanNo : mgf.scanCollection.GetMS2DescendingArray()) {
                        ScanData Scan = mgf.scanCollection.GetScan(scanNo);
                        Scan.MsLevel = 2;
                        Scan.MGFTitle = Scan.MGFTitle.split("_")[Scan.MGFTitle.split("_").length - 2];
                        Scan.Normalization();
                        Scan.Data.Finalize();
                        scanCollection.AddScan(Scan);
                    }
                } else if (files.toLowerCase().endsWith(".pkl")) {
                    PKLScanParser pkl = new PKLScanParser(folder + "/" + files);
                    pkl.scan.MsLevel = 1;
                    pkl.scan.MGFTitle = pkl.scan.MGFTitle.split("_")[pkl.scan.MGFTitle.split("_").length - 2];
                    pkl.scan.ScanNum=ScanNo;
                    //pkl.scan.Normalization();
                    pkl.scan.RetentionTime=cycletime*ScanNo;
                    ScanNo++;
                    pkl.scan.Data.Finalize();
                    SortedMap.put(pkl.scan.MGFTitle, pkl.scan);
                    //scanCollection.AddScan(pkl.scan);
                }
            }
        }
        
        for(ScanData scan : SortedMap.values()){
            scanCollection.AddScan(scan);
        }        
    }
}
