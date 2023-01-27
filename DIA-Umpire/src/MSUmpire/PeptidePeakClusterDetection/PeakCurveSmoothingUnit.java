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
package MSUmpire.PeptidePeakClusterDetection;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.PeakDataStructure.PeakCurve;
import java.util.ArrayList;
import java.util.concurrent.Callable;

/**
 * Peak shape smoothing process thread unit
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PeakCurveSmoothingUnit implements Callable<ArrayList<PeakCurve>> {

    PeakCurve curve;
    boolean export;

    InstrumentParameter parameter;

    public PeakCurveSmoothingUnit(PeakCurve curve, InstrumentParameter para) {
        this.curve = curve;
        this.parameter = para;
    }

    @Override
    public ArrayList<PeakCurve> call() {
        final ArrayList<PeakCurve> ResultCurves;
        //If we want to split multimodal peak curve by CWT
        if (parameter.DetectByCWT) {
            curve.DoBspline();
            curve.DetectPeakRegion();
            ResultCurves = curve.SeparatePeakByRegion(parameter.SNThreshold);
            curve=null;
        }
        else{
            curve.DoBspline();
            ResultCurves=new ArrayList<>();
            ResultCurves.add(curve);
        }
        /**
         * the for loop below does the work of public void ClearRawPeaks() in MSUmpire.​PeptidePeakClusterDetection.​PDHandlerBase.
         */
        for (final PeakCurve peakCurve : ResultCurves) {
            peakCurve.CalculateMzVar();
            peakCurve.StartRT();
            peakCurve.EndRT();
            peakCurve.ReleaseRawPeak();
        }

        return ResultCurves;
    }
}
