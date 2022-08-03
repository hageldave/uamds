package uamds.vis;

import java.awt.Font;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Set;
import java.util.function.DoubleBinaryOperator;
import java.util.stream.Collectors;

import org.ejml.data.DMatrixRMaj;

import hageldave.jplotter.canvas.JPlotterCanvas;
import hageldave.jplotter.misc.DefaultGlyph;
import hageldave.jplotter.renderables.Lines;
import hageldave.jplotter.renderables.Points;
import hageldave.jplotter.renderables.Text;
import hageldave.jplotter.renderers.CoordSysRenderer;
import uamds.optimization.ejml.MatCalcEJML;
import uamds.optimization.generic.numerics.MatCalc;
import uamds.other.NRV;

public class DistributionPlot<M> extends IsoLinesPlot {
	
	// percentiles: 25%-75% at 0.675 std | 10%-90% at 1.282 std | 2.5%-97.5% at 1.960 std
	static final double[] gaussStdDevs = {0.675, 1.282 , 1.960};
	static final double[] gaussPercentiles = Arrays.stream(gaussStdDevs).map(std->new NRV<DMatrixRMaj>(new MatCalcEJML(), 2).pdf().evalAt(std,0)*6.2831853).toArray();
	
	public ArrayList<Points> means = new ArrayList<>();
	public ArrayList<Points> samples = new ArrayList<>();
	public ArrayList<Lines> axes = new ArrayList<>();
	public ArrayList<Text> texts = new ArrayList<>();
	public ArrayList<String> labels = new ArrayList<>();
	protected final MatCalc<M> mc;
	protected int numRings;
	protected boolean axisEnabled=true;
	
	public DistributionPlot(MatCalc<M> mc, int numRings) {
		super();
		this.mc=mc;
		this.numRings=numRings;
	}
	
	public DistributionPlot(MatCalc<M> mc) {
		this(mc, 2);
	}
	
	public DistributionPlot<M> enableAxes(boolean enable){
		this.axisEnabled=enable;
		for(var p : means)
			p.hide(enable);
		for(var a : axes)
			a.hide(!enable);
		return this;
	}
	
	public void ensureCapacity(int n) {
		while(means.size() < n)
			addEmpty();
	}
	
	@Override
	public int addEmpty() {
		Points points = new Points(DefaultGlyph.CROSS);
		means.add(points);
		Points samplepoints = new Points(DefaultGlyph.SQUARE_F).setGlobalScaling(0.5).setGlobalAlphaMultiplier(0.5);
		samples.add(samplepoints);
		Lines lines = new Lines();
		axes.add(lines);
		completeRendererFG.addItemToRender(points);
		completeRendererFG.addItemToRender(lines);
		completeRendererFG.addItemToRender(samplepoints);
		Text txt = new Text("", 12/*16*/, Font.PLAIN);
		texts.add(txt);
		labels.add("");
		completeRendererFG.addItemToRender(txt);
		enableAxes(axisEnabled);
		return super.addEmpty();
	}
	
	public DistributionPlot(MatCalc<M> mc, JPlotterCanvas canvas, CoordSysRenderer csr) {
		super(canvas, csr);
		this.mc = mc;
	}
	
	public int addDistribution(NRV<M> nrv, int color) {
		int idx = addEmpty();
		updateDistribution(idx, nrv, color);
		return idx;
	}
	
	public void purgeDistribution(int idx) {
			means.get(idx).removeAllPoints();
			texts.get(idx).setTextString("").setOrigin(0, 0);
			purgeLines(idx);
	}
	
	public void updateDistribution(int idx, NRV<M> nrv, int color) {
		updateDistribution(idx, nrv, color, ""+idx);
	}
	
	public void updateDistribution(int idx, NRV<M> nrv, int color, String label) {
		M axes = nrv.axes();
		M bounds = mc.concatHorz(mc.addColVec(mc.scale(axes, numRings), nrv.mean), mc.addColVec(mc.scale(axes, -numRings), nrv.mean));			
//				DoubleMatrix.concatHorizontally(
//				axes.mul( 2).addColumnVector(nrv.mean), 
//				axes.mul(-2).addColumnVector(nrv.mean));
		M min = mc.rowMins(bounds); // bounds.rowMins();
		M size = mc.sub(mc.rowMaxs(bounds), min); //bounds.rowMaxs().sub(min);
		Rectangle2D boundingbox = new Rectangle2D.Double(mc.get(min,0), mc.get(min,1), mc.get(size,0), mc.get(size,1));
		
		var pdf = nrv.pdf();
		double toUnit = 1/pdf.evalAt(pdf.mu);
		M vec = mc.zeros(2);// DoubleMatrix.zeros(2);
		DoubleBinaryOperator fx = (x,y)->pdf.evalAt(mc.set_inp(mc.set_inp(vec,0,x),1,y))*toUnit;
		
//		M drawnSamples = nrv.drawSamples(1000);
		
		hageldave.jplotter.util.Utils.execOnAWTEventDispatch(()->{
			Points points = means.get(idx);
			points.removeAllPoints().addPoint(mc.get(pdf.mu,0), mc.get(pdf.mu,1))
			.setColor(color)
			.setPickColor(this.idx2pickcolor.get(idx));
			
//			Points samplePoints = samples.get(idx);
//			samplePoints.removeAllPoints();
//			for(int i=0; i<1000;i++) {
//				samplePoints.addPoint(mc.get(drawnSamples, i, 0), mc.get(drawnSamples, i, 1))
//				.setColor(color);
//			}
			
			Lines lines = this.axes.get(idx);
			lines.removeAllSegments();
			lines.addSegment(mc.get(pdf.mu,0), mc.get(pdf.mu,1), mc.get(bounds, 0, 0), mc.get(bounds, 1, 0));
			lines.addSegment(mc.get(pdf.mu,0), mc.get(pdf.mu,1), mc.get(bounds, 0, 1), mc.get(bounds, 1, 1));
			lines.addSegment(mc.get(pdf.mu,0), mc.get(pdf.mu,1), mc.get(bounds, 0, 2), mc.get(bounds, 1, 2));
			lines.addSegment(mc.get(pdf.mu,0), mc.get(pdf.mu,1), mc.get(bounds, 0, 3), mc.get(bounds, 1, 3));
			lines.getSegments().forEach(seg->seg.setColor(color).setPickColor(this.idx2pickcolor.get(idx)));

			this.labels.set(idx, label);
			Text txt = texts.get(idx);
			txt.setTextString(label);
			txt.setBackground((color&0x00ffffff)|(0x66000000));
			txt.setPickColor(this.idx2pickcolor.get(idx));
			txt.setOrigin(new Point2D.Double(mc.get(pdf.mu,0), mc.get(pdf.mu,1)));
		});
		
		// percentiles: 25%-75% at 0.675 std | 10%-90% at 1.282 std | 2.5%-97.5% at 1.960 std
		updateIsoLines(idx, fx, boundingbox, Arrays.copyOf(gaussPercentiles, numRings), new int[]{color, (color&0x00ffffff)|0x99000000, (color&0x00ffffff)|0x44000000});
	}
	
	@Override
	public synchronized void highLight(int... idx) {
		super.highLight(idx);
		for(var t : texts) {
			completeRendererBG.text.removeItemToRender(t);
			completeRendererFG.text.removeItemToRender(t);
		}
		for(var m : means) {
			completeRendererBG.points.removeItemToRender(m);
			completeRendererFG.points.removeItemToRender(m);
		}
		for(var a : axes) {
			completeRendererBG.lines.removeItemToRender(a);
			completeRendererFG.lines.removeItemToRender(a);
		}
		Set<Integer> toHighLight = Arrays.stream(idx).boxed().collect(Collectors.toSet());
		if(toHighLight.isEmpty()) {
			// reset
			for(int i=0; i<means.size(); i++) {
				Text t = texts.get(i);
				t.hide(false);
				completeRendererFG.addItemToRender(t);
				
				Points m = means.get(i);
				m.setGlobalSaturationMultiplier(1.0);
				m.setGlobalAlphaMultiplier(1.0);
				completeRendererFG.addItemToRender(m);
				
				Lines a = axes.get(i);
				a.setGlobalSaturationMultiplier(1.0);
				a.setGlobalAlphaMultiplier(1.0);
				completeRendererFG.addItemToRender(a);
			}
		} else {
			for(int i=0; i<means.size(); i++) {
				Text t = texts.get(i);
				Points m = means.get(i);
				Lines a = axes.get(i);
				if(toHighLight.contains(i)) {
					completeRendererFG.addItemToRender(t).addItemToRender(m).addItemToRender(a);
					t.hide(false);
					m.setGlobalSaturationMultiplier(1.0);
					m.setGlobalAlphaMultiplier(1.0);
					a.setGlobalSaturationMultiplier(1.0);
					a.setGlobalAlphaMultiplier(1.0);
				} else {
					completeRendererBG.addItemToRender(t).addItemToRender(m);
					t.hide(true);
					m.setGlobalSaturationMultiplier(0.2);
					m.setGlobalAlphaMultiplier(0.2);
					a.setGlobalSaturationMultiplier(0.2);
					a.setGlobalAlphaMultiplier(0.2);
				}
			}
		}
		canvas.scheduleRepaint();
	}
	
}
