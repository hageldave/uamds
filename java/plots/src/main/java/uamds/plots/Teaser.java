package uamds.plots;

import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.geom.Rectangle2D;
import java.util.Arrays;
import java.util.Random;
import java.util.function.IntUnaryOperator;
import java.util.stream.IntStream;

import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import hageldave.jplotter.canvas.BlankCanvasFallback;
import hageldave.jplotter.canvas.JPlotterCanvas;
import hageldave.jplotter.color.ColorMap;
import hageldave.jplotter.color.DefaultColorMap;
import hageldave.jplotter.color.SimpleColorMap;
import hageldave.jplotter.renderables.Lines;
import hageldave.jplotter.renderables.Text;
import hageldave.jplotter.renderers.LinesRenderer;
import hageldave.jplotter.renderers.TextRenderer;
import hageldave.jplotter.util.Utils;
import uamds.UAMDS;
import uamds.misc.Utils2;
import hageldave.optisled.ejml.MatCalcEJML;
import hageldave.optisled.generic.numerics.MatCalc;
import uamds.other.NRV;
import uamds.other.NRVSet;
import hageldave.utils.Ref;
import uamds.vis.DistributionPlot;
import uamds.vis.PerDimDistributionPlot;

public class Teaser {
	
	static int compW = 350;
	static int distrH = compW-8;
	static int dataH = compW*55/100;
	
	public static void main(String[] args) {
		MatCalc<?> mc = new MatCalcEJML();
		genTeaserFig(mc);
	}
	
	public static <M> void genTeaserFig(MatCalc<M> mc) {
		NRV<M>[] distributions;
		String[] featureLabels = null;
		IntUnaryOperator colorProvider = null;
		final int n;
		final int d;
		
		{ /* load data set */
			distributions = getTeaserDistribs(mc);
			int[] colors = DefaultColorMap.Q_9_SET1.getColors();
			int yellow = colors[5];
			colors[5] = colors[6];
			colors[6] = yellow;
			ColorMap cmap = new SimpleColorMap(colors);
			colorProvider = (i)->cmap.getColor(i%cmap.numColors());
			n = distributions.length;
			d = distributions[0].d;
			featureLabels = IntStream.range(0, d).mapToObj(i->"dim"+(i+1)).toArray(String[]::new);
		}
		// move mean of means to origin
		M center = Arrays.stream(distributions).map(rv -> rv.mean).reduce(mc::add).get();
		mc.scale_inp(center,1.0/n);
		Arrays.stream(distributions).forEach(rv -> mc.sub_inp(rv.mean,center));
		
		NRVSet<M> data0 = new NRVSet<>();
		NRVSet<M> data1 = new NRVSet<>();
		NRVSet<M> data2 = new NRVSet<>();
		for(int i=0; i<n; i++) {
			NRV<M> nrv = distributions[i];
			data0.add(nrv.scaleCov(0.0));
			data1.add(nrv.scaleCov(i==4 ? 1.0:0.0));
			data2.add(nrv);
		}
		// create initial guess for UAMDS
		M[][] init = mc.matArray(2, n);
		Random rand = new Random();
//		long seed = rand.nextLong();
//		System.out.println(Long.toHexString(seed));
//		rand.setSeed(seed);
		rand.setSeed(0x2d455607f3c50d74L);
		for(int i=0; i<n; i++) {
			init[0][i] = mc.sub(mc.rand(2, d, rand),0.5);
			init[1][i] = mc.sub(mc.rand(2, 1, rand),0.5);
		}
		// calc initial UAMDS, few iterations for zero variance, then some on the full variance
		Ref<M[][]> result = new Ref<>();
		for(int i=0; i<100; i++) {
			new UAMDS<>(mc).calculateProjection(data0, init, result, 100);
			init = result.get();
		}
		for(int i=0; i<100; i++) {
			new UAMDS<>(mc).calculateProjection(data2, init, result, 100);
			init = result.get();
		}
		// now calculate final result for data0
		// calc UAMDS for data0
		NRVSet<M> projection0=null;
		for(int i=0; i<100; i++) {
			projection0 = new UAMDS<>(mc).calculateProjection(data0, init=result.get(), result, 100);
		}
		// UAMDS for data1
		NRVSet<M> projection1=null;
		new UAMDS<>(mc).calculateProjection(data1, init, result, 100);
		for(int i=0; i<100; i++) {
			projection1 = new UAMDS<>(mc).calculateProjection(data1, result.get(), result, 100);
		}
		// UAMDS for data2
		NRVSet<M> projection2=null;
		new UAMDS<>(mc).calculateProjection(data2, init, result, 100);
		for(int i=0; i<10_0; i++) {
			projection2 = new UAMDS<>(mc).calculateProjection(data2, result.get(), result, 1000);
		}
		
		// visualization (2D projections)
		
		DistributionPlot<M> distr0 = new DistributionPlot<>(mc,1).enableAxes(false);
		DistributionPlot<M> distr1 = new DistributionPlot<>(mc,1).enableAxes(false);
		DistributionPlot<M> distr2 = new DistributionPlot<>(mc,1).enableAxes(false);
		for(int i=0; i<n; i++) {
			distr0.addEmpty();
			distr1.addEmpty();
			distr2.addEmpty();
			distr0.updateDistribution(i,projection0.get(i), colorProvider.applyAsInt(i),"");
			distr1.updateDistribution(i,projection1.get(i), colorProvider.applyAsInt(i),"");
			distr2.updateDistribution(i,projection2.get(i), colorProvider.applyAsInt(i),"");
		}
		// add points of zero variance (plot 0) to plot 1 and 2
		for(int i=0; i<n; i++) {
			distr1.addEmpty();
			distr2.addEmpty();
			int color = colorProvider.applyAsInt(i);
			color = Utils.averageColor(color, 0xffffffff);
			distr1.updateDistribution(i+n,projection0.get(i), color,"");
			distr2.updateDistribution(i+n,projection0.get(i), color,"");
		}
		
		for(DistributionPlot<?> dp:Arrays.asList(distr0,distr1,distr2)) {
			dp.coordsys.setCoordinateView(-5, -5, 5, 5);
			dp.coordsys.setPaddingBot(0).setPaddingLeft(2).setPaddingRight(2).setPaddingTop(0);
			dp.coordsys.setxAxisLabel("").setyAxisLabel("");
			dp.canvas.asComponent().setPreferredSize(new Dimension(compW, distrH));
		}
		
		// visualization (dataset whisker plots)
		
		PerDimDistributionPlot<M> pddp0 = new PerDimDistributionPlot<>(mc, d, colorProvider, featureLabels);
		PerDimDistributionPlot<M> pddp1 = new PerDimDistributionPlot<>(mc, d, colorProvider, featureLabels);
		PerDimDistributionPlot<M> pddp2 = new PerDimDistributionPlot<>(mc, d, colorProvider, featureLabels);
		pddp0.data.set(data0);
		pddp1.data.set(data1);
		pddp2.data.set(data2);
		for(PerDimDistributionPlot<?> pddp : Arrays.asList(pddp0, pddp1, pddp2)) {
			pddp.setView(new Rectangle2D.Double(0, -10, 1, 20));
			Arrays.stream(pddp.plots).forEach(plot->{
				plot.coordsys.setPaddingBot(-10).setPaddingLeft(0).setPaddingRight(-10).setPaddingTop(1);
			});
			pddp.setPreferredSize(new Dimension(compW, dataH));
		}
		
		// save visualizations to file
		
//		pddp0.display("data no variance");
//		pddp1.display("data single variance");
//		pddp2.display("data full variance");
//		distr0.display("uamds no variance");
//		distr1.display("uamds single variance");
//		distr2.display("uamds full variance");
		
		Container allPlots = new Container();
		allPlots.setBackground(Color.white);
		Container col0 = new Container();
		Container col1 = new Container();
		Container col2 = new Container();
		allPlots.setLayout(new BoxLayout(allPlots, BoxLayout.X_AXIS));
		Arrays.asList(col0,col1,col2).forEach(col->col.setLayout(new BoxLayout(col, BoxLayout.Y_AXIS)));
		col0.add(pddp0);
		col1.add(pddp1);
		col2.add(pddp2);
		col0.add(distr0.canvas.asComponent());
		col1.add(distr1.canvas.asComponent());
		col2.add(distr2.canvas.asComponent());
		allPlots.add(mkSideLabels(new Dimension(30, col0.getPreferredSize().height)));
		allPlots.add(col0);
		allPlots.add(mkSeparator(new Dimension(2, col0.getPreferredSize().height)));
		allPlots.add(col1);
		allPlots.add(mkSeparator(new Dimension(2, col0.getPreferredSize().height)));
		allPlots.add(col2);
		
		JFrame frame = Utils2.createJFrameWithBoilerPlate("teaser");
		frame.setPreferredSize(null);
		frame.setContentPane(allPlots);
		SwingUtilities.invokeLater(()->{
			frame.pack();
			frame.setVisible(true);
		});
		
		
		long time = System.currentTimeMillis();
		while(time + 500 > System.currentTimeMillis()) {
			Thread.yield(); // idle for half a second to give plots time to show and render
		}
		SwingUtilities.invokeLater(()->{
			Utils2.exportSVG(allPlots, "teaser");
			System.out.println("exported");
		});
		
	}
	
	public static Container mkSeparator(Dimension prefSize) {
		JPlotterCanvas canvas = new BlankCanvasFallback();
		LinesRenderer renderer = new LinesRenderer();
		Lines lines = new Lines();
		lines.setGlobalThicknessMultiplier(2);
		lines.addSegment(0, 1, 0, 0);
		lines.setStrokePattern(0xf0f0);
		renderer.addItemToRender(lines);
		renderer.setView(new Rectangle2D.Double(-.5, 0, 1, 1));
		canvas.setRenderer(renderer);
		Container c = new Container();
		c.setLayout(new BoxLayout(c, BoxLayout.Y_AXIS));
		c.add(canvas.asComponent());
		c.setPreferredSize(prefSize);
		return c;
	}
	
	public static Container mkSideLabels(Dimension prefSize) {
		JPlotterCanvas canvas = new BlankCanvasFallback();
		TextRenderer renderer = new TextRenderer();
		Text lblDat = new Text("Dataset", 28, Font.PLAIN);
		lblDat.setAngle(Math.PI/2).setOrigin(29, distrH+(dataH-(int)lblDat.getBounds().getWidth())/2);
		
		Text lblProj = new Text("UAMDS", 28, Font.PLAIN);
		lblProj.setAngle(Math.PI/2).setOrigin(29, (distrH-(int)lblProj.getBounds().getWidth())/2);
		
		renderer.addItemToRender(lblDat);
		renderer.addItemToRender(lblProj);
		canvas.setRenderer(renderer);
		Container c = new Container();
		c.setLayout(new BoxLayout(c, BoxLayout.Y_AXIS));
		c.add(canvas.asComponent());
		c.setPreferredSize(prefSize);
		return c;
	}
	
	@SuppressWarnings("unchecked")
	static <M> NRV<M>[] getTeaserDistribs(MatCalc<M> mc){
		return new NRV[] {
				new NRV<>(mc, mc.vecOf(0.615, -1.852, 0.028, -0.430), mc.matOf(4, 5.099, 0.000, 0.000, 0.000, 0.000, 0.039, 0.000, 0.000, 0.000, 0.000, 1.241, 0.000, 0.000, 0.000, 0.000, 5.428)),
				new NRV<>(mc, mc.vecOf(0.934, -1.360, 1.678, 0.264), mc.matOf(4, 1.111, 0.000, 0.000, 0.000, 0.000, 0.132, 0.000, 0.000, 0.000, 0.000, 0.021, 0.000, 0.000, 0.000, 0.000, 4.565)),
				new NRV<>(mc, mc.vecOf(-2.340, 2.034, -0.666*1.3, 2.584*0.5), mc.matOf(4, 5.764, 0.000, 0.000, 0.000, 0.000, 2.782, 0.000, 0.000, 0.000, 0.000, 8.146, 0.000, 0.000, 0.000, 0.000, 8.752)),
				new NRV<>(mc, mc.vecOf(-0.774, -1.583, -2.382, -1.007), mc.matOf(4, 1.626, 0.000, 0.000, 0.000, 0.000, 0.400, 0.000, 0.000, 0.000, 0.000, 1.328, 0.000, 0.000, 0.000, 0.000, 0.140)),
				new NRV<>(mc, mc.vecOf(-1.645, -0.326, -0.769, -1.027), mc.matOf(4, 8.778, 0.000, 0.000, 0.000, 0.000, 2.035, 0.000, 0.000, 0.000, 0.000, 0.337, 0.000, 0.000, 0.000, 0.000, 3.387)),
				new NRV<>(mc, mc.vecOf(0.730, 2.257, -0.656, -2.485), mc.matOf(4, 4.636, 0.000, 0.000, 0.000, 0.000, 0.534, 0.000, 0.000, 0.000, 0.000, 1.021, 0.000, 0.000, 0.000, 0.000, 3.034)),
		};
	}
	
}