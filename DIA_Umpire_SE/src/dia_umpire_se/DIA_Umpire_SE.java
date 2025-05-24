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
package dia_umpire_se;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.BaseDataStructure.UmpireInfo;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.DIA.DIAPack;
import MSUmpire.Utility.ConsoleLogger;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Locale;
import java.util.concurrent.TimeUnit;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class DIA_Umpire_SE {

    /**
     * @param args the command line arguments DIA_Umpire parameterfile
     */
    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);

        Thread.setDefaultUncaughtExceptionHandler((t, e) -> {
            e.printStackTrace();
            System.exit(1);
        });

        if (org.burningwave.core.assembler.StaticComponentContainer.Modules != null)
            org.burningwave.core.assembler.StaticComponentContainer.Modules.exportAllToAll(); // for FST serialization to work an Java 16 and above
        System.out.println("=================================================================================================");
        System.out.println("DIA-Umpire signal extraction analysis  (version: " + UmpireInfo.GetInstance().Version + ")");
        System.out.println("(c) University of Michigan");
        if (args.length < 2 || args.length > 3) {
            System.out.println("command format error, the correct format is: java -jar -Xmx8G DIA_Umpire_SE.jar mzMXL_file diaumpire_se.params");
            System.out.println("To fix DIA setting, use : java -jar -Xmx8G DIA_Umpire_SE.jar mzMXL_file diaumpire_se.params -f");
            return;
        }
        {
            //Define logger level for console
            ConsoleLogger.SetConsoleLogger(Level.INFO);
            //Define logger level and file path for text log file
            ConsoleLogger.SetFileLogger(Level.DEBUG, FilenameUtils.getFullPath(args[0]) + "diaumpire_se.log");
        }

        boolean Fix = false;
        boolean Resume=false;
        
        if (args.length == 3 && args[2].equals("-f")) {
            Fix = true;
        }
        String parameterfile = args[1];
        String MSFilePath = args[0];
        LogManager.getRootLogger().info("Version: " + UmpireInfo.GetInstance().Version);
        LogManager.getRootLogger().info("Parameter file:" + parameterfile);
        LogManager.getRootLogger().info("Spectra file:" + MSFilePath);
        BufferedReader reader = new BufferedReader(new FileReader(parameterfile));

        String line = "";
        InstrumentParameter param = new InstrumentParameter(InstrumentParameter.InstrumentType.TOF5600);
        param.DetermineBGByID = false;
        param.EstimateBG = true;
        int NoCPUs = 2;

        SpectralDataType.DataType dataType = SpectralDataType.DataType.DIA_F_Window;
        String WindowType = "";
        int WindowSize = 25;

        ArrayList<XYData> WindowList = new ArrayList<>();

        boolean ExportPrecursorPeak = false;
        boolean ExportFragmentPeak = false;

        //<editor-fold defaultstate="collapsed" desc="Read parameter file">
        while ((line = reader.readLine()) != null) {
            LogManager.getRootLogger().info(line);
            if (!"".equals(line) && !line.startsWith("#")) {
                //System.out.println(line);
                if (line.equals("==window setting begin")) {
                    while (!(line = reader.readLine()).equals("==window setting end")) {
                        if (!"".equals(line)) {
                            WindowList.add(new XYData(Float.parseFloat(line.split("\t")[0]), Float.parseFloat(line.split("\t")[1])));
                        }
                    }
                    continue;
                }
                if (line.split("=").length < 2) {
                    continue;
                }
                String type = line.split("=")[0].trim();
                if (type.startsWith("para.")) {
                    type = type.replace("para.", "SE.");
                }
                String value = line.split("=")[1].trim();
                switch (type) {
                    case "Thread": {
                        NoCPUs = Integer.parseInt(value);
                        break;
                    }
                    case "ExportPrecursorPeak": {
                        ExportPrecursorPeak = Boolean.parseBoolean(value);
                        break;
                    }
                    case "ExportFragmentPeak": {
                        ExportFragmentPeak = Boolean.parseBoolean(value);
                        break;
                    }

                    //<editor-fold defaultstate="collapsed" desc="instrument parameters">
                    case "RPmax": {
                        param.PrecursorRank = Integer.parseInt(value);
                        break;
                    }
                    case "RFmax": {
                        param.FragmentRank = Integer.parseInt(value);
                        break;
                    }
                    case "CorrThreshold": {
                        param.CorrThreshold = Float.parseFloat(value);
                        break;
                    }
                    case "DeltaApex": {
                        param.ApexDelta = Float.parseFloat(value);
                        break;
                    }
                    case "RTOverlap": {
                        param.RTOverlapThreshold = Float.parseFloat(value);
                        break;
                    }
                    case "BoostComplementaryIon": {
                        param.BoostComplementaryIon = Boolean.parseBoolean(value);
                        break;
                    }
                    case "AdjustFragIntensity": {
                        param.AdjustFragIntensity = Boolean.parseBoolean(value);
                        break;
                    }
                    case "Q1":
                        param.Q1 = Boolean.parseBoolean(value);
                        break;
                    case "Q2":
                        param.Q2 = Boolean.parseBoolean(value);
                        break;
                    case "Q3":
                        param.Q3 = Boolean.parseBoolean(value);
                        break;
                    case "SE.MS1PPM": {
                        param.MS1PPM = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MS2PPM": {
                        param.MS2PPM = Float.parseFloat(value);
                        break;
                    }
                    case "SE.SN": {
                        param.SNThreshold = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MS2SN": {
                        param.MS2SNThreshold = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MinMSIntensity": {
                        param.MinMSIntensity = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MinMSMSIntensity": {
                        param.MinMSMSIntensity = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MinRTRange": {
                        param.MinRTRange = Float.parseFloat(value);
                        break;
                    }
                    case "SE.MaxNoPeakCluster": {
                        param.MaxNoPeakCluster = Integer.parseInt(value);
                        param.MaxMS2NoPeakCluster = Integer.parseInt(value);
                        break;
                    }
                    case "SE.MinNoPeakCluster": {
                        param.MinNoPeakCluster = Integer.parseInt(value);
                        param.MinMS2NoPeakCluster = Integer.parseInt(value);
                        break;
                    }
                    case "SE.MinMS2NoPeakCluster": {
                        param.MinMS2NoPeakCluster = Integer.parseInt(value);
                        break;
                    }
                    case "SE.MaxCurveRTRange": {
                        param.MaxCurveRTRange = Float.parseFloat(value);
                        break;
                    }
                    case "SE.Resolution": {
                        param.Resolution = Integer.parseInt(value);
                        break;
                    }
                    case "SE.RTtol": {
                        param.RTtol = Float.parseFloat(value);
                        break;
                    }
                    case "SE.NoPeakPerMin": {
                        param.NoPeakPerMin = Integer.parseInt(value);
                        break;
                    }
                    case "SE.StartCharge": {
                        param.StartCharge = Integer.parseInt(value);
                        break;
                    }
                    case "SE.EndCharge": {
                        param.EndCharge = Integer.parseInt(value);
                        break;
                    }
                    case "SE.MS2StartCharge": {
                        param.MS2StartCharge = Integer.parseInt(value);
                        break;
                    }
                    case "SE.MS2EndCharge": {
                        param.MS2EndCharge = Integer.parseInt(value);
                        break;
                    }
                    case "SE.NoMissedScan": {
                        param.NoMissedScan = Integer.parseInt(value);
                        break;
                    }
                    case "SE.Denoise": {
                        param.Denoise = Boolean.valueOf(value);
                        break;
                    }
                    case "SE.EstimateBG": {
                        param.EstimateBG = Boolean.valueOf(value);
                        break;
                    }
                    case "SE.RemoveGroupedPeaks": {
                        param.RemoveGroupedPeaks = Boolean.valueOf(value);
                        break;
                    }
                    case "SE.MinFrag": {
                        param.MinFrag = Integer.parseInt(value);
                        break;
                    }
                    case "SE.IsoPattern": {
                        param.IsoPattern = Float.valueOf(value);
                        break;
                    }
                    case "SE.StartRT": {
                        param.startRT = Float.valueOf(value);
                        break;
                    }
                    case "SE.EndRT": {
                        param.endRT = Float.valueOf(value);
                        break;
                    }
                    case "SE.RemoveGroupedPeaksRTOverlap": {
                        param.RemoveGroupedPeaksRTOverlap = Float.valueOf(value);
                        break;
                    }
                    case "SE.RemoveGroupedPeaksCorr": {
                        param.RemoveGroupedPeaksCorr = Float.valueOf(value);
                        break;
                    }
                    case "SE.MinMZ": {
                        param.MinMZ = Float.valueOf(value);
                        break;
                    }
                    case "SE.MinPrecursorMass":{
                        param.MinPrecursorMass = Float.valueOf(value);
                        break;
                    }
                    case "SE.MaxPrecursorMass":{
                        param.MaxPrecursorMass = Float.valueOf(value);
                        break;
                    }
                    case "SE.IsoCorrThreshold": {
                        param.IsoCorrThreshold = Float.valueOf(value);
                        break;
                    }
                    case "SE.MassDefectFilter": {
                        param.MassDefectFilter = Boolean.parseBoolean(value);
                        break;
                    }
                    case "SE.MassDefectOffset": {
                        param.MassDefectOffset = Float.valueOf(value);
                        break;
                    }
                    

//</editor-fold>//</editor-fold>
                    
                    case "WindowType": {
                        WindowType = value;
                        switch (WindowType) {
                            case "SWATH": {
                                dataType = SpectralDataType.DataType.DIA_F_Window;
                                break;
                            }
                            case "V_SWATH": {
                                dataType = SpectralDataType.DataType.DIA_V_Window;
                                break;
                            }
                            case "MSX": {
                                dataType = SpectralDataType.DataType.MSX;
                                break;
                            }
                            case "MSE": {
                                dataType = SpectralDataType.DataType.MSe;
                                break;
                            }
                        }
                        break;
                    }
                    case "WindowSize": {
                        WindowSize = Integer.parseInt(value);
                        break;
                    }
                }
            }
        }
//</editor-fold>

        try {
            File MSFile = new File(MSFilePath);
            if (MSFile.exists()) {
                long time = System.currentTimeMillis();
                LogManager.getRootLogger().info("=================================================================================================");
                LogManager.getRootLogger().info("Processing " + MSFilePath + "....");
                
                //Initialize a DIA file data structure                
                DIAPack DiaFile = new DIAPack(MSFile.getAbsolutePath(), NoCPUs);
                DiaFile.Resume=Resume;
                DiaFile.SetDataType(dataType);
                DiaFile.SetParameter(param);
                
                //Set DIA isolation window setting
                if (dataType == SpectralDataType.DataType.DIA_F_Window) {
                    DiaFile.SetWindowSize(WindowSize);
                } else if (dataType == SpectralDataType.DataType.DIA_V_Window) {
                    for (XYData window : WindowList) {
                        DiaFile.AddVariableWindow(window);
                    }
                }
                DiaFile.SaveDIASetting();                
                DiaFile.SaveParams();
                
                if (Fix) {
                    DiaFile.FixScanidx();
                    return;
                }
                DiaFile.ExportPrecursorPeak = ExportPrecursorPeak;
                DiaFile.ExportFragmentPeak = ExportFragmentPeak;
                LogManager.getRootLogger().info("Module A: Signal extraction");
                //Start DIA signal extraction process to generate pseudo MS/MS files
                DiaFile.process();
                time = System.currentTimeMillis() - time;
                LogManager.getRootLogger().info(MSFilePath + " processed time:" + String.format("%d hour, %d min, %d sec", TimeUnit.MILLISECONDS.toHours(time), TimeUnit.MILLISECONDS.toMinutes(time) - TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(time)), TimeUnit.MILLISECONDS.toSeconds(time) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(time))));
            }else{
                throw new RuntimeException("file: “"+MSFile+"” does not exist!");
            }
            LogManager.getRootLogger().info("Job complete");
            LogManager.getRootLogger().info("=================================================================================================");

        } catch (Exception e) {
            LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(e));
            throw e;
        }
    }
}
