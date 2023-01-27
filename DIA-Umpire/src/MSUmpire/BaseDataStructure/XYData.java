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

import java.io.Serializable;

/*
 * Two dimensional data
 */
/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class XYData implements Comparable<XYData>, Serializable {

    private static final long serialVersionUID = 973492749274921L;

    private float X, Y;

//    private synchronized void writeObject(java.io.ObjectOutputStream stream) throws java.io.IOException {        
//        stream.defaultWriteObject();
//        stream.writeFloat(xydata[0]);        
//        stream.writeFloat(xydata[1]);        
//    }
//    private void readObject(java.io.ObjectInputStream in) throws ClassNotFoundException, IOException {
//        in.defaultReadObject();
//        xydata[0]=in.readFloat();
//        xydata[1]=in.readFloat();
//    }
    public XYData(float x, float y) {
        setX(x);
        setY(y);
    }

    @Override
    public int compareTo(XYData o) {
        return Float.compare(o.getX(), getX());
    }

    /**
     * @return the X
     */
    public float getX() {
        return X;
    }

    /**
     * @param X the X to set
     */
    public void setX(float X) {
        this.X = X;
    }

    /**
     * @return the Y
     */
    public float getY() {
        return Y;
    }

    /**
     * @param Y the Y to set
     */
    public void setY(float Y) {
        this.Y = Y;
    }

    public XYData cloneXYData() {
        return new XYData(getX(), getY());
    }
}
