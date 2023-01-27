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
import MSUmpire.SpectralProcessingModule.BackgroundDetector;
import MSUmpire.SpectralProcessingModule.ScanPeakGroup;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

/**
 * MGF file parser
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class MGFParser {

    public ScanCollection scanCollection = null;
    public String filename;
    public InstrumentParameter parameter;
    public int TotalScan;
    public int NoCPUs = 4;
    public SpectralDataType.DataType datatype;

    public MGFParser(String filename, InstrumentParameter parameter, SpectralDataType.DataType datatype, int NoCPUs) throws IOException {
        this.filename = filename;
        scanCollection = new ScanCollection(parameter.Resolution);
        scanCollection.Filename = filename;
        this.parameter = parameter;
        this.datatype = datatype;
        this.NoCPUs = NoCPUs;
        //Parse();
    }

    public MGFParser(String filename) throws IOException {
        parameter = new InstrumentParameter(InstrumentParameter.InstrumentType.Q_TOF);
        this.filename = filename;
        scanCollection = new ScanCollection(parameter.Resolution);
        scanCollection.Filename = filename;
        this.datatype = SpectralDataType.DataType.DDA;
        this.NoCPUs = 12;
        //Parse();
    }

    public void GetAllScanCollectionDDA() throws FileNotFoundException, IOException {
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        String line = "";
        String[] Values = null;
        int ScanNum = 1;

        while ((line = reader.readLine()) != null) {
            if (line.trim().startsWith("BEGIN IONS")) {
                ScanData scan = new ScanData();
                while (!(line = reader.readLine()).trim().startsWith("END IONS")) {
                    if (line.trim().startsWith("PEPMASS=")) {
                        scan.PrecursorMz = Float.parseFloat(line.trim().subSequence(8, line.trim().length()).toString());
                    }
                    if (line.trim().startsWith("CHARGE=")) {
                        scan.PrecursorCharge = Integer.parseInt(line.trim().subSequence(7, line.trim().length() - 1).toString());
                    }
                    if (line.trim().startsWith("RTINSECONDS=")) {
                        scan.RetentionTime = Float.parseFloat(line.trim().subSequence(12, line.trim().length()).toString()) / 60f;
                    }                    
                    if (line.trim().startsWith("TITLE=")) {
                        scan.MGFTitle = line.trim().subSequence(6, line.trim().length()).toString().replace(",", "_");
                        if (scan.MGFTitle.contains(" RT:")) {
                            scan.RetentionTime = Float.parseFloat(scan.MGFTitle.substring(scan.MGFTitle.indexOf(" RT:") + 4).split(" ")[0]);
                        }
                        if (scan.MGFTitle.contains("=Scan:")) {
                            scan.ScanNum = Integer.parseInt(scan.MGFTitle.substring(scan.MGFTitle.indexOf("=Scan:") + 6).split(" ")[0]);
                            ScanNum = scan.ScanNum;
                        }
                    }
                    if ((Values = line.split(" ")).length == 2) {
                        scan.AddPoint(Float.parseFloat(Values[0]), Float.parseFloat(Values[1]));
                    } else if ((Values = line.split("\t")).length == 2) {
                        scan.AddPoint(Float.parseFloat(Values[0]), Float.parseFloat(Values[1]));
                    }
                }
                if (scan.ScanNum == 0) {
                    scan.ScanNum = ScanNum++;
                }
                if (scan.PrecursorCharge == 0) {
                    scan.MsLevel = 1;
                } else {
                    scan.MsLevel = 2;
                }
                scan.centroided = false;
                scan.Data.Finalize();
                scan.background = 0f;
                if (parameter.EstimateBG) {
                    BackgroundDetector detector=new BackgroundDetector(scan);
                    detector.DetermineConstantBackground();
                } else {
                    if (scan.MsLevel == 1) {
                        scan.background = parameter.MinMSIntensity;
                    }
                    if (scan.MsLevel == 2) {
                        scan.background = parameter.MinMSMSIntensity;
                    }
                }

                if (parameter.Denoise) {
                    scan.RemoveSignalBelowBG();
                }

                if (!scan.centroided) {
                    scan.Centroiding(parameter.Resolution, scan.background);
                }
                if (parameter.Deisotoping && scan.MsLevel == 1) {
                   ScanPeakGroup scanpeak= new ScanPeakGroup(scan, parameter);
                   scanpeak.Deisotoping();
                }
                scanCollection.AddScan(scan);
            }
        }
        reader.close();
    }

    private void Parse() throws FileNotFoundException, IOException {
        //        BEGIN IONS
//PEPMASS=820.998855732003
//CHARGE=1+
//RTINSECONDS=200
//TITLE=Elution from: 0.14 to 0.14   period: 0   experiment: 2 cycles:  1
//200.9942 2.3857
//354.9856 2.3857
//370.9314 5.1571
//388.9714 9.6857
//390.9608 2.7429
//END IONS
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        String line = "";
        String[] Values = null;
        int ScanNum = 1;

        while ((line = reader.readLine()) != null) {
            if (line.trim().startsWith("BEGIN IONS")) {
                ScanData scan = new ScanData();
                scan.ScanNum = ScanNum++;
                scan.MsLevel = 2;
                while (!(line = reader.readLine()).trim().startsWith("END IONS")) {
                    if (line.trim().startsWith("PEPMASS=")) {
                        scan.PrecursorMz = Float.parseFloat(line.trim().subSequence(8, line.trim().length()).toString());
                    }
                    if (line.trim().startsWith("CHARGE=")) {
                        scan.PrecursorCharge = Integer.parseInt(line.trim().subSequence(7, line.trim().length() - 1).toString());
                    }
                    if (line.trim().startsWith("RTINSECONDS=")) {
                        scan.RetentionTime = Float.parseFloat(line.trim().subSequence(12, line.trim().length()).toString());
                    }
                    if (line.trim().startsWith("TITLE=")) {
                        scan.MGFTitle = line.trim().subSequence(6, line.trim().length()).toString().replace(",", "_");
                    }
                    if ((Values = line.split(" ")).length == 2) {
                        scan.AddPoint(Float.parseFloat(Values[0]), Float.parseFloat(Values[1]));
                    } else if ((Values = line.split("\t")).length == 2) {
                        scan.AddPoint(Float.parseFloat(Values[0]), Float.parseFloat(Values[1]));
                    }
                }
                scan.Data.Finalize();
                scanCollection.AddScan(scan);
            }
        }
        reader.close();
    }
}
