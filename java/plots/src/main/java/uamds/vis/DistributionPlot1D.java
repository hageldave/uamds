package uamds.vis;

import java.awt.Dimension;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.IntUnaryOperator;
import java.util.stream.IntStream;

import org.ejml.data.DMatrixRMaj;

import uamds.datasets.StudentGrades;
import uamds.misc.Utils2;
import hageldave.jplotter.color.DefaultColorMap;
import hageldave.jplotter.color.SimpleColorMap;
import hageldave.jplotter.renderables.Lines;
import hageldave.jplotter.renderables.Text;
import hageldave.jplotter.renderers.CompleteRenderer;
import hageldave.optisled.ejml.MatCalcEJML;
import hageldave.optisled.generic.numerics.MatCalc;
import uamds.other.Utils;
import uamds.other.MultivariateGaussian;
import uamds.other.NRV;
import uamds.other.NRVSet;
import hageldave.utils.Ref;

public class DistributionPlot1D<M> extends CoordSysDisplay {

	MatCalc<M> mc;
	public Ref<NRVSet<M>> distribs1D = new Ref<>();
	public String[] labels = null;
	ArrayList<Lines> lines = new ArrayList<>();
	ArrayList<Text> texts = new ArrayList<>();
	CompleteRenderer content = new CompleteRenderer();
	IntUnaryOperator colorForDistrib = (i)->0xff222222;
	public boolean horizontal;
	public boolean normalize;
	public boolean boxplot;
	
	public DistributionPlot1D(MatCalc<M> mc, boolean opengl, boolean horizontal) {
		super(opengl);
		this.mc = mc;
		this.coordsys.setContent(content);
		this.horizontal = horizontal;
		distribs1D.addListener(this::ondistribChange);
		canvas.asComponent().setPreferredSize(new Dimension(400, 120));
		Utils2.addExportMenuOnMiddleClick(canvas);
		
		this.coordsys.addCoordinateViewListener((src, view) -> {
			if(horizontal && view.getY() != 0) {
				coordsys.setCoordinateView(new Rectangle2D.Double(view.getMinX(), 0, view.getWidth(), view.getHeight()));
			}
			if(!horizontal && view.getX() != 0) {
				coordsys.setCoordinateView(new Rectangle2D.Double(0, view.getMinY(), view.getWidth(), view.getHeight()));
			}
		});
	}
	
	public Rectangle2D getContentBounds() {
		if(distribs1D.isNull() || distribs1D.get().isEmpty()) {
			return new Rectangle2D.Double(0,0,1,1);
		} else {
			Rectangle2D bounds = lines.get(0).getBounds();
			for(Lines l:lines) {
				bounds = l.getBounds().createUnion(bounds);
			}
			return bounds;
		}
	}
	
	void ondistribChange(NRVSet<M> distribs) {
		while(lines.size() < distribs.size()) {
			Lines l = new Lines();
			lines.add(l);
			content.addItemToRender(l);
			Text txt = new Text("", 12, Font.PLAIN);
			texts.add(txt);
			content.addItemToRender(txt);
		}
		
		int numSamples = 64;
		double toUnit = 1.0/(numSamples-1);
		double maxVal = distribs.stream()
				.filter(nrv->mc.get(nrv.cov,0) > 1e-6)
				.mapToDouble(nrv->{
				return nrv.pdf().evalAt(nrv.mean);
		}).max().orElse(1.0);
		
		for(int i=0; i<distribs.size(); i++) {
			NRV<M> nrv = distribs.get(i);
			double std = Math.sqrt(mc.get(nrv.cov, 0));
			double mu = mc.get(nrv.mean, 0);
			int color = colorForDistrib.applyAsInt(i);
			if(boxplot) {
//				lines.get(i).setVertexRoundingEnabled(true);
				// percentiles: 25%-75% at 0.675 std | 10%-90% at 1.282 std | 2.5%-97.5% at 1.960 std
				double q1 = mu-0.675*std;
				double q3 = mu+0.675*std;
				double iqr = q3-q1;
				double whisk1 = q1-1.5*iqr;
				double whisk2 = q3+1.5*iqr;
				
				double c = i+0.5;
				double e = 0.4;
				if(normalize) {
					c /= distribs.size();
					e /= distribs.size();
				}
				lines.get(i).removeAllSegments();
				if(horizontal) {
					// mean
					lines.get(i).addSegment(mu, c-e, mu, c+e);
					if(std > 1e-6) {
						// box
						lines.get(i).addSegment(q1, c-e, q1, c+e);
						lines.get(i).addSegment(q3, c-e, q3, c+e);
						lines.get(i).addSegment(q1, c-e, q3, c-e);
						lines.get(i).addSegment(q1, c+e, q3, c+e);
						// whisker 1
						lines.get(i).addSegment(q1, c, whisk1, c);
						lines.get(i).addSegment(whisk1, c-e, whisk1, c+e);
						// whisker 2
						lines.get(i).addSegment(q3, c, whisk2, c);
						lines.get(i).addSegment(whisk2, c-e, whisk2, c+e);
					}
				} else {
					// mean
					lines.get(i).addSegment(c-e, mu, c+e, mu);
					if(std > 1e-6) {
						// box
						lines.get(i).addSegment(c-e, q1, c+e, q1);
						lines.get(i).addSegment(c-e, q3, c+e, q3);
						lines.get(i).addSegment(c-e, q3, c-e, q1);
						lines.get(i).addSegment(c+e, q3, c+e, q1);
						// whisker 1
						lines.get(i).addSegment(c, whisk1, c, q1);
						lines.get(i).addSegment(c-e, whisk1, c+e, whisk1);
						// whisker 2
						lines.get(i).addSegment(c, whisk2, c, q3);
						lines.get(i).addSegment(c-e, whisk2, c+e, whisk2);
					}
				}
				// color
				lines.get(i).getSegments().forEach(seg->seg.setColor(color));
			} else { /* density plot */
//				lines.get(i).setVertexRoundingEnabled(false);
				double[] x,y;
				if(std > 1e-6) {
					MultivariateGaussian<M> pdf = nrv.pdf();
					x = IntStream.range(0, numSamples).mapToDouble(j->(mu-3.5*std)+(j*toUnit)*(7*std)).toArray();
					y = Arrays.stream(x).map(pdf::evalAt).toArray();
					if(normalize) {
						for(int j=0; j<y.length; j++)
							y[j] *= (1/maxVal)*0.95;
					}
				} else {
					x = new double[]{mu,mu};
					y = new double[]{0,1};
				}
				if(horizontal) {
					lines.get(i).removeAllSegments().addLineStrip(x, y).forEach(seg->seg.setColor(color));
				} else {
					lines.get(i).removeAllSegments().addLineStrip(y, x).forEach(seg->seg.setColor(color));
				}
				if(labels != null && labels.length > i) {
					Text txt = texts.get(i);
					txt.setTextString(labels[i]);
					if(horizontal) {
						txt.setOrigin(new Point2D.Double(mu, Arrays.stream(y).max().getAsDouble()*1.1));
					} else {
						txt.setOrigin(new Point2D.Double(Arrays.stream(y).max().getAsDouble()*1.1, mu));
					}
				}
			}
		}
		
		canvas.scheduleRepaint();
	}
	
	public void autoSetView() {
		Rectangle2D bounds = getContentBounds();
		if(horizontal) {
			coordsys.setCoordinateView(bounds.getX(), bounds.getY(), bounds.getWidth(), bounds.getHeight()*1.1);
		} else {
			coordsys.setCoordinateView(bounds.getX(), bounds.getY(), bounds.getWidth()*1.2, bounds.getHeight());
		}
	}

	public static void main(String[] args) {
		MatCalc<DMatrixRMaj> mc = new MatCalcEJML();
		DistributionPlot1D<DMatrixRMaj> plot = new DistributionPlot1D<>(mc, false, true);
		plot.coordsys.setCoordinateView(0, 0, 20, 1);
		NRVSet<DMatrixRMaj> data = new NRVSet<>();
		
		double[][] mD = StudentGrades.mark2distrib_2;
		
		data.add( new NRV<DMatrixRMaj>(mc, mc.vecOf(mD[0][0]), mc.vecOf(mD[0][1])));
		data.add( new NRV<DMatrixRMaj>(mc, mc.vecOf(mD[1][0]), mc.vecOf(mD[1][1])));
		data.add( new NRV<DMatrixRMaj>(mc, mc.vecOf(mD[2][0]), mc.vecOf(mD[2][1])));
		data.add( new NRV<DMatrixRMaj>(mc, mc.vecOf(mD[3][0]), mc.vecOf(mD[3][1])));
		data.add( new NRV<DMatrixRMaj>(mc, mc.vecOf(mD[4][0]), mc.vecOf(mD[4][1])));
		data.add( new NRV<DMatrixRMaj>(mc, mc.vecOf(mD[5][0]), mc.vecOf(mD[5][1])));
		
		SimpleColorMap coolwarm = DefaultColorMap.D_COOL_WARM.reversed();
		SimpleColorMap cmap = new SimpleColorMap(
				coolwarm.interpolate(0),
				coolwarm.interpolate(0.15),
				coolwarm.interpolate(0.3),
				coolwarm.interpolate(0.7),
				coolwarm.interpolate(0.85),
				coolwarm.interpolate(1)
				);
		plot.colorForDistrib = cmap::getColor;
		plot.coordsys.setxAxisLabel("grade").setyAxisLabel("density");
		plot.labels = "VB,B,FB,FG,G,VG".split(",");
		plot.distribs1D.set(data);
		plot.setupAndShow();
	}
	
}
