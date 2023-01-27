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
package MSUmpire.PeakDataStructure;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.XYData;
import java.io.Serializable;
import java.util.ArrayList;

/*
 *This class implements the Continuous Wavelet Transform (CWT), Mexican Hat,
 * over raw datapoints of a certain spectrum. After get the spectrum in the
 * wavelet's time domain, we use the local maxima to detect possible peaks in
 * the original raw datapoints.
 * Described in Tautenhahn, R., Bottcher, C. & Neumann, S. 
 * Highly sensitive feature detection for high resolution LC/MS. 
 * BMC Bioinformatics 9, 504 (2008).
 */
public class WaveletMassDetector implements Serializable{

    /**
     * Parameters of the wavelet, NPOINTS is the number of wavelet values to use
     * The WAVELET_ESL & WAVELET_ESL indicates the Effective Support boundaries
     */
    private final double NPOINTS;
    private final int WAVELET_ESL = -5;
    private final int WAVELET_ESR = 5;
    static boolean waveletDebug = false;
    private final InstrumentParameter parameter;
//    public ArrayList<XYData> DataPoint;
    public final float[] DataPoint;
    double waveletWindow = 0.3;
    private final double[] MEXHAT;
    double NPOINTS_half;

//    public WaveletMassDetector(InstrumentParameter parameter, ArrayList<XYData> DataPoint, int NoPoints) {
    public WaveletMassDetector(InstrumentParameter parameter, final float[] DataPoint, int NoPoints) {
        this.parameter = parameter;
        this.DataPoint = DataPoint;
        this.NPOINTS = NoPoints;

        double wstep = ((WAVELET_ESR - WAVELET_ESL) / NPOINTS);
        MEXHAT = new double[(int) NPOINTS];

        double waveletIndex = WAVELET_ESL;
        for (int j = 0; j < NPOINTS; j++) {
            // Pre calculate the values of the wavelet
            MEXHAT[j] = cwtMEXHATreal(waveletIndex, waveletWindow, 0.0);
            waveletIndex += wstep;
        }

        NPOINTS_half = NPOINTS / 2;
        d = (int) NPOINTS / (WAVELET_ESR - WAVELET_ESL);
    }
    int d;
    //public ArrayList<XYData>[] waveletCWT;
    
    //List of peak ridge (local maxima)
    public ArrayList<XYData>[] PeakRidge;

    public void Run() {

        //"Intensities less than this value are interpreted as noise",                
        //"Scale level",
        //"Number of wavelet'scale (coeficients) to use in m/z peak detection"
        //"Wavelet window size (%)",
        //"Size in % of wavelet window to apply in m/z peak detection");        
//        int maxscale = (int) (Math.max(Math.min((DataPoint.get(DataPoint.size() - 1).getX() - DataPoint.get(0).getX()), parameter.MaxCurveRTRange), 0.5f) * parameter.NoPeakPerMin / (WAVELET_ESR + WAVELET_ESR));
        int maxscale = (int) (Math.max(Math.min((DataPoint[2*(DataPoint.length/2 - 1)] - DataPoint[0]), parameter.MaxCurveRTRange), 0.5f) * parameter.NoPeakPerMin / (WAVELET_ESR + WAVELET_ESR));

        //waveletCWT = new ArrayList[15];
        PeakRidge = new ArrayList[maxscale];
        //XYData maxint = new XYData(0f, 0f);
        for (int scaleLevel = 0; scaleLevel < maxscale; scaleLevel++) {
//            ArrayList<XYData> wavelet = performCWT(scaleLevel * 2 + 5);
            final float[] wavelet = performCWT(scaleLevel * 2 + 5);
            PeakRidge[scaleLevel] = new ArrayList<>();
            //waveletCWT[scaleLevel] = wavelet;
//            XYData lastpt = wavelet.get(0);
            int lastptidx = 0;
//            XYData localmax = null;
            int localmaxidx = -1;
//            XYData startpt = wavelet.get(0);
            int startptidx = 0;

            boolean increasing = false;
            boolean decreasing = false;
//            XYData localmaxint = null;
            int localmaxintidx = -1;

//            for (int cwtidx = 1; cwtidx < wavelet.size(); cwtidx++) {
            for (int cwtidx = 1; cwtidx < wavelet.length/2; cwtidx++) {
//                XYData CurrentPoint = wavelet.get(cwtidx);
                final float CurrentPointY = wavelet[2*cwtidx+1],
                        lastptY=wavelet[2*lastptidx+1],
                        startptY=wavelet[2*startptidx+1],
                        localmaxY=localmaxidx==-1?Float.NaN:wavelet[2*localmaxidx+1];
//                if (CurrentPoint.getY() > lastpt.getY()) {//the peak is increasing
                if (CurrentPointY> lastptY) {//the peak is increasing
                    if (decreasing) {//first increasing point, last point was a possible local minimum
                        //check if the peak was symetric
//                        if (localmax != null && (lastpt.getY() <= startpt.getY() || Math.abs(lastpt.getY() - startpt.getY()) / localmax.getY() < parameter.SymThreshold)) {
                        if (localmaxidx != -1 && (lastptY<= startptY|| Math.abs(lastptY- startptY) / localmaxY < parameter.SymThreshold)) {
//                            PeakRidge[scaleLevel].add(localmax);
                            PeakRidge[scaleLevel].add(new XYData(wavelet[2*localmaxidx],wavelet[2*localmaxidx+1]));
//                            localmax = CurrentPoint;
                            localmaxidx = cwtidx;
//                            startpt = lastpt;
                            startptidx = lastptidx;
                        }
                    }
                    increasing = true;
                    decreasing = false;
//                } else if (CurrentPoint.getY() < lastpt.getY()) {//peak decreasing
                } else if (CurrentPointY < lastptY) {//peak decreasing
                    if (increasing) {//first point decreasing, last point was a possible local maximum
//                        if (localmax == null || localmax.getY() < lastpt.getY()) {
                        if (localmaxidx == -1 || localmaxY < lastptY) {
//                            localmax = lastpt;
                            localmaxidx = lastptidx;
                        }
                    }
                    decreasing = true;
                    increasing = false;
                }
//                lastpt = CurrentPoint;
                lastptidx = cwtidx;
                final float localmaxintY=localmaxintidx == -1 ? Float.NaN : wavelet[2*localmaxintidx+1];
//                if (localmaxint == null || CurrentPoint.getY() > localmaxint.getY()) {
                if (localmaxintidx == -1 || CurrentPointY > localmaxintY) {
//                    localmaxint = CurrentPoint;
                    localmaxintidx = cwtidx;
                }
//                if (cwtidx == wavelet.size() - 1 && decreasing) {
                if (cwtidx == wavelet.length/2 - 1 && decreasing) {
//                    if (localmax != null && (CurrentPoint.getY() <= startpt.getY() || Math.abs(CurrentPoint.getY() - startpt.getY()) / localmax.getY() < parameter.SymThreshold)) {
//                    final float startptY=wavelet[2*startptidx+1];
                    if (localmaxidx != -1 && (CurrentPointY <= startptY || Math.abs(CurrentPointY - startptY) / localmaxY < parameter.SymThreshold)) {
                        final XYData localmax = new XYData(wavelet[2*localmaxidx], wavelet[2*localmaxidx+1]);
                        PeakRidge[scaleLevel].add(localmax);
                    }
                }
            }

            if (!waveletDebug) {
//                wavelet.clear();
                //wavelet = null;
            }
        }
    }

    /**
     * Perform the CWT over raw data points in the selected scale level
     *
     *
     */
//    private ArrayList<XYData> performCWT(int scaleLevel) {
    private float[] performCWT(int scaleLevel) {
//        int length = DataPoint.size();
        int length = DataPoint.length/2;
//        ArrayList<XYData> cwtDataPoints = new ArrayList<XYData>();
        final float[] cwtDataPoints = new float[length*2];

        int a_esl = scaleLevel * WAVELET_ESL;
        int a_esr = scaleLevel * WAVELET_ESR;
        double sqrtScaleLevel = Math.sqrt(scaleLevel);
        for (int dx = 0; dx < length; dx++) {
            /*
             * Compute wavelet boundaries
             */
            int t1 = a_esl + dx;
            if (t1 < 0) {
                t1 = 0;
            }
            int t2 = a_esr + dx;
            if (t2 >= length) {
                t2 = (length - 1);
            }

            /*
             * Perform convolution
             */
            float intensity = 0f;
            for (int i = t1; i <= t2; i++) {
                int ind = (int) (NPOINTS_half) + (d * (i - dx) / scaleLevel);
                if (ind < 0) {
                    ind = 0;
                }
                if (ind >= NPOINTS) {
                    ind = (int) NPOINTS - 1;
                }
//                if(i<0 || ind<0){
//                    System.out.print("");
//                }
//                intensity += DataPoint.get(i).getY() * MEXHAT[ind];
                intensity += DataPoint[2*i+1] * MEXHAT[ind];
            }
            intensity /= sqrtScaleLevel;
            // Eliminate the negative part of the wavelet map
            if (intensity < 0) {
                intensity = 0;
            }
//            cwtDataPoints.add(new XYData(DataPoint.get(dx).getX(), intensity));
            cwtDataPoints[2*dx]= DataPoint[2*dx];
            cwtDataPoints[2*dx+1]= intensity;
        }
        return cwtDataPoints;
    }

    /**
     * This function calculates the wavelets's coefficients in Time domain
     *
     * @param x Step of the wavelet
     * @param window Window Width of the wavelet
     * @param b Offset from the center of the peak
     */
    private double cwtMEXHATreal(double x, double window, double b) {
        /*
         * c = 2 / ( sqrt(3) * pi^(1/4) )
         */
        double c = 0.8673250705840776;
        double TINY = 1E-200;
        double x2;

        if (window == 0.0) {
            window = TINY;
        }
        //x-b=t
        //window=delta
        x = (x - b) / window;
        x2 = x * x;
        return c * (1.0 - x2) * Math.exp(-x2 / 2);
    }
    /**
     * This function searches for maximums from wavelet data points
     */
}
