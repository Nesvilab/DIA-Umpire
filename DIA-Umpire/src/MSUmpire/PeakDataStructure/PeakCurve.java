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
import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.BaseDataStructure.XYZData;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.PriorityQueue;
import org.eclipse.collections.impl.list.mutable.primitive.FloatArrayList;


/**
 * Single m/z trace peak curve
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PeakCurve implements Serializable  {
    public static class XYZDataList{
        public XYZDataList(){}
        final FloatArrayList xyzdata=new FloatArrayList();
        public void addXYZ(XYZData d){
            xyzdata.ensureCapacity(xyzdata.size()+3);
            xyzdata.add(d.getX());
            xyzdata.add(d.getY());
            xyzdata.add(d.getZ());
        }
        public void add(float x,float y, float z){
            xyzdata.ensureCapacity(xyzdata.size()+3);
            xyzdata.add(x);
            xyzdata.add(y);
            xyzdata.add(z);
        }
        public XYZData get(int i){return new XYZData(getXat(i),getYat(i),getZat(i));}
        public float getXat(int i){return xyzdata.get(3*i);}
        public float getYat(int i){return xyzdata.get(3*i+1);}
        public float getZat(int i){return xyzdata.get(3*i+2);}

        public void setXat(int i,float f){xyzdata.set(3*i,f);}
        public void setYat(int i,float f){xyzdata.set(3*i+1,f);}
        public void setZat(int i,float f){xyzdata.set(3*i+2,f);}

        public int size(){return xyzdata.size()/3;}
        public void clear(){xyzdata.clear();}
    }
    private static final long serialVersionUID = 6498163564821L;

//    private ArrayList<XYZData> PeakList;
    private  XYZDataList PeakList;
    //X: retention time
    //Y: m/z
    //Z: intensity
    private XYPointCollection SmoothData;
    //X: retention time
    //Y: intensity
    //XIC
    private float startint = 0f;
    public int Index;
    private float endrt = -1f;
    private float startrt = -1f;
    public int StartScan=-1;
    public int EndScan=-1;
    private float TotalIntMzF;
    private float TotalIntF;
    public float TargetMz;
    public float ApexInt;
    public float minIntF = Float.POSITIVE_INFINITY;
    public float ApexRT;
    public float MaxCorr = 0f;
    public boolean CheckState = false;
    public float ConflictCorr = 0f;
    public boolean Grouped = false;
//    public transient HashSet<Integer> ChargeGrouped=new HashSet<>();
    public transient short ChargeGrouped=0;
    public float MzVar = -1f;
    public transient SortedRidgeCollectionClass PeakRidgeList;
    public transient WaveletMassDetector waveletMassDetector;
//    private transient ArrayList<XYZData> PeakRegionList;
    private transient XYZDataList PeakRegionList;
//    public transient ArrayList<Float> RegionRidge;
    public transient FloatArrayList RegionRidge;
//    private transient ArrayList<ArrayList<Float>> NoRidgeRegion;
    private transient ArrayList<FloatArrayList> NoRidgeRegion;
    public InstrumentParameter parameter;

    //using B-spline to generate smoothed peak signals
    public void DoBspline() {
//        for (XYZData point : PeakList) {
        for (int i=0; i<PeakList.size();++i) {
//            XYData pt = new XYData(point.getX(), point.getZ());
            XYData pt = new XYData(PeakList.getXat(i), PeakList.getZat(i));
            SmoothData.AddPoint(pt);
        }
        SmoothData.Data.Finalize();// to sorted array
        SmoothData = new Bspline().Run(SmoothData, (int) Math.max((RTWidth() * parameter.NoPeakPerMin), PeakList.size()), 2);
    }
    public void DoBspline_() {
        final long[] unslong = new long[PeakList.size()];
        for (int i=0; i<PeakList.size();++i)
            unslong[i]=Bspline.ff_to_long(PeakList.getXat(i), PeakList.getZat(i));

        Arrays.sort(unslong);
        final float[] sorted=new float[unslong.length*2];
        for(int i=0;i<sorted.length/2;++i){
            sorted[2*i]=Bspline.get_X(unslong[i]);
            sorted[2*i+1]=Bspline.get_Y(unslong[i]);
        }
        final float[] smoothed=Bspline.Run__(sorted, (int) Math.max((RTWidth() * parameter.NoPeakPerMin), PeakList.size()), 2);
        for(int i=0;i<smoothed.length/2;++i)
            this.SmoothData.AddPoint(new XYData(smoothed[2*i],smoothed[2*i+1]));
        this.SmoothData.Data.Finalize();
    }
    public void DoInterpolation() {
//        for (XYZData point : PeakList) {
        for (int i=0; i<PeakList.size();++i) {
//            XYData pt = new XYData(point.getX(), point.getZ());
            XYData pt = new XYData(PeakList.getXat(i), PeakList.getZat(i));
            SmoothData.AddPoint(pt);
        }
        LinearInterpolation interpo = new LinearInterpolation();
        SmoothData = interpo.Run(SmoothData, (int) Math.max((RTWidth() * parameter.NoPeakPerMin), PeakList.size()));
        interpo = null;
    }

    public void AddConflictScore(float corr) {
        synchronized (this) {
            ConflictCorr += corr;
        }
    }

    public float GetRawSNR() {
        return ApexInt / minIntF;
    }

    public void SetRTs(float StartRT, float EndRT) {
        startrt = StartRT;
        endrt = EndRT;
    }

    //Detect peak region using CWT based on smoothed peak signals
    public void DetectPeakRegion() {
        PeakRidgeList = new SortedRidgeCollectionClass();
//        PeakRegionList = new ArrayList<>();
        PeakRegionList = new XYZDataList();
        NoRidgeRegion = new ArrayList<>();
        if (RTWidth() * parameter.NoPeakPerMin < 1) {
            return;
        }

//        ArrayList<XYData> PeakArrayList = new ArrayList<>();
        final float[] PeakArrayList = new float[SmoothData.PointCount()*2];
        for (int i = 0; i < SmoothData.PointCount(); i++) {
//            PeakArrayList.add(new XYData(SmoothData.Data.get(i).getX(), SmoothData.Data.get(i).getY()));
            PeakArrayList[2*i]=SmoothData.Data.get(i).getX();
            PeakArrayList[2*i+1]=SmoothData.Data.get(i).getY();
        }
        //Start CWT process
        waveletMassDetector = new WaveletMassDetector(parameter, PeakArrayList, (int) (RTWidth() * parameter.NoPeakPerMin));
        waveletMassDetector.Run();

        int maxScale = waveletMassDetector.PeakRidge.length - 1;

        float[] DisMatrixF = null;//new float[PeakRidgeList.size()*waveletMassDetector.PeakRidge[maxScale].size()*2];
        //trace peak ridge from maximum wavelet scale to minimum scale
        for (int i = maxScale; i >= 0; i--) {
            //Get peak ridge list (maximum RT points given a CWT scale
            ArrayList<XYData> PeakRidgeArray = waveletMassDetector.PeakRidge[i];

            if (PeakRidgeArray == null) {
                maxScale = i;
                continue;
            }
            if (PeakRidgeArray.isEmpty()) {
                continue;
            }

            //RT distance matrix between the groupped peak riges and peak ridges extracted from current CWT scale
            final int r=PeakRidgeList.size(), c=PeakRidgeArray.size();
//            float[][] DisMatrixF = new float[PeakRidgeList.size()][PeakRidgeArray.size()];
//            final float[] DisMatrixF = new float[r*c];
            if(DisMatrixF==null || r*c>DisMatrixF.length)
                DisMatrixF = new float[r*c*2];

            for (int k = 0; k < PeakRidgeList.size(); k++) {///For each existing peak ridge line
                for (int l = 0; l < PeakRidgeArray.size(); l++) {
//                    DisMatrixF[k][l] = Math.abs(PeakRidgeList.get(k).RT - PeakRidgeArray.get(l).getX());
                    DisMatrixF[k*c+l] = Math.abs(PeakRidgeList.get(k).RT - PeakRidgeArray.get(l).getX());
                }
            }

            boolean conti = true;
            ArrayList<XYData> RemovedRidgeList = new ArrayList<>();
            while (conti) {
                float closest = Float.MAX_VALUE;
                int ExistingRideIdx = -1;
                int PeakRidgeInx = -1;
                for (int k = 0; k < PeakRidgeList.size(); k++) {
                    for (int l = 0; l < PeakRidgeArray.size(); l++) {
                        {
//                            if (DisMatrixF[k][l] < closest) {
                            if (DisMatrixF[k*c+l] < closest) {
//                                closest = DisMatrixF[k][l];
                                closest = DisMatrixF[k*c+l];
                                ExistingRideIdx = k;
                                PeakRidgeInx = l;
                            }
                        }
                    }
                }

                if (closest < Float.MAX_VALUE && closest <= parameter.MinRTRange) {
                    PeakRidge ridge = PeakRidgeList.remove(ExistingRideIdx);
                    ridge.lowScale = i;
                    ridge.ContinuousLevel++;
                    XYData nearestRidge = PeakRidgeArray.get(PeakRidgeInx);
                    ridge.RT = nearestRidge.getX();
                    PeakRidgeList.add(ridge);
                    RemovedRidgeList.add(nearestRidge);
                    for (int k = 0; k < PeakRidgeList.size(); k++) {
                        DisMatrixF[k*c+PeakRidgeInx] = Float.MAX_VALUE;
//                        DisMatrixF[k][PeakRidgeInx] = Float.MAX_VALUE;
                    }
                    for (int l = 0; l < PeakRidgeArray.size(); l++) {
                        DisMatrixF[ExistingRideIdx*c+l] = Float.MAX_VALUE;
//                        DisMatrixF[ExistingRideIdx][l] = Float.MAX_VALUE;
                    }
                } else {
                    conti = false;
                }
            }

            PeakRidgeArray.removeAll(RemovedRidgeList);

            RemovedRidgeList.clear();
            RemovedRidgeList = null;
            ArrayList<PeakRidge> removelist = new ArrayList<>();
            for (int k = 0; k < PeakRidgeList.size(); k++) {
                PeakRidge existridge = PeakRidgeList.get(k);
                if (existridge.lowScale - i > 2 && existridge.ContinuousLevel < maxScale / 2) {
                    removelist.add(existridge);
                }
            }
            for (int k = 0; k < removelist.size(); k++) {
                PeakRidgeList.remove(removelist.get(k));
            }
            removelist.clear();
            removelist = null;
            if (i > maxScale / 2) {
                for (XYData ridge : PeakRidgeArray) {
                    PeakRidge newRidge = new PeakRidge();
                    newRidge.RT = ridge.getX();
                    newRidge.lowScale = i;
                    newRidge.ContinuousLevel++;
                    newRidge.intensity = SmoothData.GetPoinByXCloset(newRidge.RT).getY();
                    PeakRidgeList.add(newRidge);
                }
            }
            PeakRidgeArray.clear();
            PeakRidgeArray = null;
        }

        if (PeakRidgeList.size() <= 1) {
            PeakRegionList.add(SmoothData.Data.get(0).getX(), ApexRT, SmoothData.Data.get(SmoothData.PointCount() - 1).getX());
//            ArrayList<Float> RidgeRTs = new ArrayList<>();
            final FloatArrayList RidgeRTs = new FloatArrayList();
            RidgeRTs.add(ApexRT);
            NoRidgeRegion.add(RidgeRTs);
        }
        if (PeakRidgeList.size() > 1) {
            XYData[] ValleyPoints = new XYData[PeakRidgeList.size() + 1];
            ValleyPoints[0] = SmoothData.Data.get(0).cloneXYData();
            PeakRidge currentridge = PeakRidgeList.get(0);
            XYData localmin = new XYData(-1f, Float.MAX_VALUE);
            int startidx = SmoothData.GetLowerIndexOfX(currentridge.RT);

            for (int j = 1; j < PeakRidgeList.size(); j++) {
                PeakRidge nextridge = PeakRidgeList.get(j);
                for (int i = startidx; i < SmoothData.Data.size(); i++) {
                    XYData point = SmoothData.Data.get(i);
                    if (point.getX() > currentridge.RT && point.getX() < nextridge.RT) {
                        if (localmin.getY() > point.getY()) {
                            localmin = point.cloneXYData();
                        }
                    }
                    if (point.getX() >= nextridge.RT) {
                        startidx = i;
                        break;
                    }
                }
                ValleyPoints[j] = localmin;
                localmin = new XYData(-1f, Float.MAX_VALUE);
                currentridge = nextridge;
            }
            ValleyPoints[PeakRidgeList.size()] = SmoothData.Data.get(SmoothData.PointCount() - 1).cloneXYData();

            //Correct ridge rt and intensity
            startidx = 0;
            for (int i = 0; i < PeakRidgeList.size(); i++) {
                PeakRidge ridge = PeakRidgeList.get(i);
                for (int j = startidx; j < SmoothData.Data.size(); j++) {
                    XYData point = SmoothData.Data.get(j);
                    if (point.getX() < ValleyPoints[i + 1].getX()) {
                        if (ridge.intensity < point.getY()) {
                            ridge.intensity = point.getY();
                            ridge.RT = point.getX();
                        }
                    } else {
                        startidx = j;
                        break;
                    }
                }
            }

            //Find split points to generate peak regions
            boolean[] Splitpoints = new boolean[PeakRidgeList.size() - 1];
            int left = 0;
            int right = PeakRidgeList.size() - 1;
            FindSplitPoint(left, right, ValleyPoints, Splitpoints);

//            ArrayList<Float> RidgeRTs = new ArrayList<>();
            FloatArrayList RidgeRTs = new FloatArrayList();
            startidx = 0;
            PeakRidge maxridge = PeakRidgeList.get(0);

            for (int i = 0; i < PeakRidgeList.size() - 1; i++) {
                RidgeRTs.add(PeakRidgeList.get(i).RT);
                if (PeakRidgeList.get(i).intensity > maxridge.intensity) {
                    maxridge = PeakRidgeList.get(i);
                }
                if (Splitpoints[i]) {
                    PeakRegionList.add(ValleyPoints[startidx].getX(), maxridge.RT, ValleyPoints[i + 1].getX());
                    NoRidgeRegion.add(RidgeRTs);

                    maxridge = PeakRidgeList.get(i + 1);
//                    RidgeRTs = new ArrayList<>();
                    RidgeRTs = new FloatArrayList();
                    startidx = i + 1;
                }
            }
            RidgeRTs.add(PeakRidgeList.get(PeakRidgeList.size() - 1).RT);
            if (PeakRidgeList.get(PeakRidgeList.size() - 1).intensity > maxridge.intensity) {
                maxridge = PeakRidgeList.get(PeakRidgeList.size() - 1);
            }
            PeakRegionList.add(ValleyPoints[startidx].getX(), maxridge.RT, ValleyPoints[PeakRidgeList.size()].getX());
            RidgeRTs.trimToSize();
            NoRidgeRegion.add(RidgeRTs);
        }
        waveletMassDetector = null;
        PeakRidgeList.clear();
        PeakRidgeList = null;
    }

    private void FindSplitPoint(int left, int right, XYData[] ValleyPoints, boolean[] splitpoints) {
        for (int i = left; i < right; i++) {
            if (ValidSplitPoint(left, right, i, ValleyPoints)) {
                splitpoints[i] = true;
                FindSplitPoint(left, i, ValleyPoints, splitpoints);
                FindSplitPoint(i + 1, right, ValleyPoints, splitpoints);
                break;
            }
        }
    }

    private boolean ValidSplitPoint(int left, int right, int cut, XYData[] ValleyPoints) {

        PeakRidge leftridge = PeakRidgeList.get(left);
        PeakRidge rightridge = PeakRidgeList.get(cut + 1);

        for (int i = left; i <= cut; i++) {
            if (PeakRidgeList.get(i).intensity > leftridge.intensity) {
                leftridge = PeakRidgeList.get(i);
            }
        }
        for (int i = cut + 1; i <= right; i++) {
            if (PeakRidgeList.get(i).intensity > rightridge.intensity) {
                rightridge = PeakRidgeList.get(i);
            }
        }
        return (Math.abs(ValleyPoints[left].getY() - ValleyPoints[cut + 1].getY()) / leftridge.intensity < parameter.SymThreshold && Math.abs(ValleyPoints[cut + 1].getY() - ValleyPoints[right + 1].getY()) / rightridge.intensity < parameter.SymThreshold);
    }

    //Split if multiple peak curves are detected
    public ArrayList<PeakCurve> SeparatePeakByRegion(float SN) {

        ArrayList<PeakCurve> tempArrayList = new ArrayList<>();
        ArrayList<PeakCurve> returnArrayList = new ArrayList<>();

        //Generate a peak curve for each detected region
        for (int i = 0; i < GetPeakRegionList().size(); i++) {
            PeakCurve peakCurve = new PeakCurve(parameter);
            peakCurve.RegionRidge = NoRidgeRegion.get(i);
            tempArrayList.add(peakCurve);
//            XYZData region = GetPeakRegionList().get(i);
            float x= GetPeakRegionList().getXat(i);
            float y= GetPeakRegionList().getYat(i);
            float z= GetPeakRegionList().getZat(i);
            if (z - x > parameter.MaxCurveRTRange) {
                int leftidx = GetSmoothedList().GetLowerIndexOfX(x);
                int rightidx = GetSmoothedList().GetHigherIndexOfX(z);
                XYData left = GetSmoothedList().Data.get(leftidx);
                XYData right = GetSmoothedList().Data.get(rightidx);
                while ((right.getX() - left.getX()) > parameter.MaxCurveRTRange) {
                    if (right.getX() - y <= parameter.MaxCurveRTRange / 4f) {
                        leftidx++;
                    } else if (y - left.getX() <= parameter.MaxCurveRTRange / 4f) {
                        rightidx--;
                    } else if (left.getY() < right.getY()) {
                        leftidx++;
                    } else {
                        rightidx--;
                    }
                    left = GetSmoothedList().Data.get(leftidx);
                    right = GetSmoothedList().Data.get(rightidx);
                }
//                region.setX(left.getX());
                GetPeakRegionList().setXat(i,left.getX());
//                region.setZ(right.getX());
                GetPeakRegionList().setZat(i,right.getX());
            }
        }

        //Add corresponding raw peaks
        for (int i = 0; i < GetPeakList().size(); i++) {
//            XYZData peak = GetPeakList().get(i);
            float x = GetPeakList().getXat(i);
            float y = GetPeakList().getYat(i);
            float z = GetPeakList().getZat(i);
            for (int j = 0; j < GetPeakRegionList().size(); j++) {
//                XYZData region = GetPeakRegionList().get(j);
                float regionx = GetPeakRegionList().getXat(j);
                float regionz = GetPeakRegionList().getZat(j);
                if (x >= regionx && x <= regionz) {
                    tempArrayList.get(j).AddPeak(x,y,z);
                    break;
                }
            }
        }

        //Add corresponding smoothed peaks
        for (int i = 0; i < GetSmoothedList().Data.size(); i++) {
            XYData peak = GetSmoothedList().Data.get(i);
            for (int j = 0; j < GetPeakRegionList().size(); j++) {
//                XYZData region = GetPeakRegionList().get(j);
                float regionx = GetPeakRegionList().getXat(j);
                float regionz = GetPeakRegionList().getZat(j);
                if (peak.getX() >= regionx && peak.getX() <= regionz) {
                    tempArrayList.get(j).GetSmoothedList().Data.add(peak);
                    break;
                }
            }
        }

        for (PeakCurve peak : tempArrayList) {
            if (peak.PeakList.size() > 2) {
                peak.GetSmoothedList().Data.Finalize();                
                returnArrayList.add(peak);
            }
        }        
        return returnArrayList;
    }

    public PeakCurve(InstrumentParameter parameter) {
        this.parameter = parameter;
        SmoothData = new XYPointCollection();
//        PeakList = new ArrayList<>();
        PeakList = new XYZDataList();
//        PeakRegionList = new ArrayList<>();
        PeakRegionList = new XYZDataList();
    }

    public float StartInt() {
        if (startint == 0f) {
//            startint = PeakList.get(0).getZ();
            startint = PeakList.getZat(0);
        }
        return startint;
    }

    public float StartRT() {
        if (startrt == -1f) {
            if (SmoothData != null && SmoothData.Data.size() > 0) {
                startrt = SmoothData.Data.get(0).getX();
            } else {
//                startrt = PeakList.get(1).getX();
                startrt = PeakList.getXat(1);
            }
        }
        return startrt;
    }
    float _snr = -1f;

    public float GetSNR() {
        if (_snr == -1) {
            _snr = ApexInt;
        }
        return _snr;
    }

    public float GetMaxIntensityByRegionRange(float StartRT, float EndRT) {
        float max = 0f;
        for (int j = 0; j < GetSmoothedList().PointCount(); j++) {
            XYData pt = GetSmoothedList().Data.get(j);
            if (pt.getX() >= StartRT && pt.getX() <= EndRT && pt.getY() > max) {
                max = pt.getY();
            }
        }
        return max;
    }

    private void CalculateBaseLine() {
        _baseLine = 0f;
        PriorityQueue<Float> IntensityQueue = new PriorityQueue<>();
        for (XYData point : SmoothData.Data) {
            IntensityQueue.add(point.getY());
        }

        if (IntensityQueue.size() > 10) {
            for (int i = 0; i < IntensityQueue.size() / 10; i++) {
                _baseLine += IntensityQueue.poll();
            }
            _baseLine /= (IntensityQueue.size() / 10);
        } else {
            _baseLine = IntensityQueue.poll();
        }
    }
    float _baseLine = -1f;

    public float GetBaseLine() {
        if (_baseLine == -1) {
            CalculateBaseLine();
            if (_baseLine == 0) {
                _baseLine = 1f;
            }
        }
        return _baseLine;
    }
    float _noiseLevel = -1f;

    public float GetNoiseLevel() {
        if (_noiseLevel == -1) {
            CalculateBaseLine();
        }
        return _noiseLevel;
    }

    public float EndRT() {
        if (endrt == -1f) {
            endrt = PeakList.getXat(PeakList.size() - 2);
        }
        return endrt;
    }

    public float LastScanRT() {
//        return PeakList.get(PeakList.size() - 1).getX();
        return PeakList.getXat(PeakList.size() - 1);
    }

    public XYPointCollection GetPeakCollection() {
        XYPointCollection PtCollection = new XYPointCollection();

        for (int i = 0; i < SmoothData.Data.size(); i++) {
            PtCollection.AddPoint(SmoothData.Data.get(i).getX(), SmoothData.Data.get(i).getY());
        }
        return PtCollection;
    }

    public XYPointCollection GetSmoothPeakCollection(float startRT, float endRT) {
        XYPointCollection PtCollection = new XYPointCollection();

        for (int i = 0; i < SmoothData.PointCount(); i++) {
            XYData pt = SmoothData.Data.get(i);
            if (pt.getX() > endRT) {
                break;
            } else if (pt.getX() >= startRT && pt.getX() <= endRT) {
                PtCollection.AddPoint(pt.getX(), pt.getY());
            }
        }
        return PtCollection;
    }

    public float DetermineIntByRTRange(float StartRT, float EndRT) {
        float Intensity = 0f;
        for (int j = 0; j < GetSmoothedList().PointCount(); j++) {
            XYData pt = GetSmoothedList().Data.get(j);
            if (pt.getX() >= StartRT && pt.getX() <= EndRT) {
                if (pt.getY() > Intensity) {
                    Intensity = pt.getY();
                }
            }
        }
        return Intensity;
    }

    public float RTWidth() {

        float Width = 0f;
        if (PeakList.size() > 0) {
//            Width = PeakList.get(PeakList.size() - 1).getX() - PeakList.get(0).getX();
            Width = PeakList.getXat(PeakList.size() - 1) - PeakList.getXat(0);
        } else if (SmoothData.PointCount() > 0) {
            Width = SmoothData.Data.get(SmoothData.PointCount() - 1).getX() - SmoothData.Data.get(0).getX();
        }
        return Width;
    }
//    public ArrayList<XYZData> GetPeakList() {
    public XYZDataList GetPeakList() {
        return PeakList;
    }
    public void nullifyPeakList() {
        this.PeakList=null;
    }

    public XYPointCollection GetSmoothedList() {
        return SmoothData;
    }

    public XYZDataList GetPeakRegionList() {
        return PeakRegionList;
    }

    public void ReleasePeakData() {

        this.PeakList = null;
        this.SmoothData = null;
        this.PeakRegionList = null;
        this.PeakRidgeList = null;
//        this.waveletMassDetector.DataPoint.clear();
        this.waveletMassDetector = null;
    }

    public void ReleaseRawPeak() {
        this.PeakList = null;
        this.PeakRegionList = null;
        this.PeakRidgeList = null;
        this.waveletMassDetector = null;
    }

//    public void AddPeak(XYZData xYZPoint) {
//
//        PeakList.add(xYZPoint);
//        TotalIntMzF += xYZPoint.getY() * xYZPoint.getZ() * xYZPoint.getZ();
//        TotalIntF += xYZPoint.getZ() * xYZPoint.getZ();
//        if (xYZPoint.getZ() > ApexInt) {
//            ApexInt = xYZPoint.getZ();
//            ApexRT = xYZPoint.getX();
//        }
//        if (xYZPoint.getZ() < minIntF) {
//            minIntF = xYZPoint.getZ();
//        }
//        TargetMz = TotalIntMzF / TotalIntF;
//    }
    public void AddPeak(float x, float y, float z) {

        PeakList.add(x,y,z);
        TotalIntMzF += y * z * z;
        TotalIntF += z*z;
        if (z > ApexInt) {
            ApexInt = z;
            ApexRT = x;
        }
        if (z < minIntF) {
            minIntF = z;
        }
        TargetMz = TotalIntMzF / TotalIntF;
    }

    public void CalculateMzVar() {
        MzVar = 0f;
        for (int j = 0; j < PeakList.size(); j++) {
//            MzVar += (PeakList.get(j).getX() - TargetMz) * (PeakList.get(j).getX() - TargetMz);
            MzVar += (PeakList.getXat(j) - TargetMz) * (PeakList.getXat(j) - TargetMz);
        }
        MzVar /= PeakList.size();
    }

}
