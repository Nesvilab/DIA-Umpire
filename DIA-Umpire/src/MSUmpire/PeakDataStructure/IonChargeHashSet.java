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

/**
 *
 * @author
 */
public final class IonChargeHashSet {
    private short state=0;
    public void add(final int charge){
        assert charge>0 && charge<16 : "charge must be in (0,16):\t" + charge;
        this.state |= 1 << charge;
    }
    public boolean contains(final int charge){
        assert charge>0 && charge<16 : "charge must be in (0,16):\t" + charge;
//        return (this.state & 1 << charge) != 0;
        return (this.state >>> charge & 1) != 0;
    }

    public static short add(final short state,final int charge){
        assert charge>0 && charge<16 : "charge must be in (0,16):\t" + charge;
        return (short) (state | 1 << charge);
    }
    public static boolean contains(final short state, final int charge){
        assert charge>0 && charge<16 : "charge must be in (0,16):\t" + charge;
        return (state >>> charge & 1) != 0;
    }

}
