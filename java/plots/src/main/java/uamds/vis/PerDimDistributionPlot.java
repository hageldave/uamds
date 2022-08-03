package uamds.vis;

import java.awt.Container;
import java.awt.Dimension;
import java.awt.geom.Rectangle2D;
import java.util.Arrays;
import java.util.function.IntUnaryOperator;

import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import hageldave.jplotter.coordsys.TickMarkGenerator;
import hageldave.jplotter.renderers.CoordSysRenderer;
import hageldave.jplotter.util.Pair;
import uamds.misc.Utils2;
import uamds.optimization.generic.numerics.MatCalc;
import uamds.other.NRV;
import uamds.other.NRVSet;
import uamds.other.Ref;

public class PerDimDistributionPlot<M> {

	MatCalc<M> mc;
	int d;
	public Ref<NRVSet<M>> data = new Ref<NRVSet<M>>(new NRVSet<>());
	DistributionPlot1D<M>[] plots; 
	IntUnaryOperator colorForDistrib = (i)->0xff222222;
	String[] featureLabels;

	public PerDimDistributionPlot(MatCalc<M> mc, int d, IntUnaryOperator colors, String[] featureLabels) {
		this.mc = mc;
		this.d = d;
		this.colorForDistrib = colors == null ? colorForDistrib:colors;
		this.featureLabels = featureLabels;
		setup();
		data.addListener(this::onDataChange);
	}
	
	void setup() {
		this.plots = new DistributionPlot1D[d];
		for(int i=0; i<d; i++) {
			this.plots[i] = new DistributionPlot1D<>(mc, false, false);
			plots[i].boxplot=true;
			plots[i].normalize=true;
			plots[i].colorForDistrib = this.colorForDistrib;
			plots[i].normalize=true;
			plots[i].lines.forEach(l->l.setGlobalThicknessMultiplier(2));
			CoordSysRenderer csys = plots[i].coordsys;
			TickMarkGenerator delegate = csys.getTickMarkGenerator();
			csys.setTickMarkGenerator(new TickMarkGenerator() {
				@Override
				public Pair<double[], String[]> genTicksAndLabels(double min, double max, int desiredNumTicks,
						boolean verticalAxis) {
					Pair<double[], String[]> tickslabels = delegate.genTicksAndLabels(min, max, desiredNumTicks, verticalAxis);
					if(verticalAxis)
						return tickslabels;
					double[] ticks = new double[]{-10};
					String[] labels = new String[ticks.length];
					Arrays.fill(labels, "");
					return Pair.of(ticks, labels);
				}
			});
			csys.setxAxisLabel(featureLabels == null || featureLabels.length == 0 ? "" : featureLabels[i]).setyAxisLabel("");
			csys.setPaddingRight(-10).setPaddingLeft(0);
		}
	}
	
	void onDataChange(NRVSet<M> nrvs) {
		for(int dim=0; dim<d; dim++) {
			NRVSet<M> set1D = new NRVSet<>();
			M proj = mc.set_inp(mc.zeros(1, d), 0,dim, 1.0);
			M zero = mc.zeros(1);
			for(NRV<M> nrvHi:nrvs) {
				NRV<M> nrvLo = nrvHi.affineTransform(proj, zero);
				set1D.add(nrvLo);
			}
			plots[dim].colorForDistrib = this.colorForDistrib;
			plots[dim].distribs1D.set(set1D);
//			plots[dim].lines.forEach(l->l.setGlobalThicknessMultiplier(2));
			plots[dim].coordsys.setxAxisLabel(featureLabels == null || featureLabels.length == 0 ? "" : featureLabels[dim]);
		}
	}
	
	public void setView(Rectangle2D rect) {
		for(DistributionPlot1D<M> plot: plots) {
			plot.coordsys.setCoordinateView(rect);
		}
	}
	
	
	public JFrame display(String title) {
		JFrame frame = Utils2.createJFrameWithBoilerPlate(title);
		Container container = new  Container();
		container.setLayout(new BoxLayout(container, BoxLayout.X_AXIS));
		frame.getContentPane().add(container);
		for(DistributionPlot1D<M> plot: plots) {
			container.add(plot.canvas.asComponent());
			plot.canvas.asComponent().setPreferredSize(new Dimension(400/4, 300));
		}
		SwingUtilities.invokeLater(()->{
			frame.pack();
			frame.setVisible(true);
		});
		
		return frame;
	}
	
}
