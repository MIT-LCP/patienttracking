

package com.ikarosilva.pttrack;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicLong;


public class  DistanceMatrix implements Callable<Double>{
	private static final int MAX_THREADS=Runtime.getRuntime().availableProcessors();
	private final int N,M;		//Total number of records
	private final BlockingQueue<Integer> tasks;
	private final HashMap<Integer,Integer> index;
	private AtomicLong results;  //Keep track of processed records
	private final double[][] database;

	public DistanceMatrix(double[][] database) throws InterruptedException{
		//Set task queue 
		N= database.length;
		M=database[0].length;
		tasks=new ArrayBlockingQueue<Integer>(N);
		this.database=database;
		index=new HashMap<Integer,Integer>();
		Integer ind=0;
		for(int n=0;n<N;n++){
			tasks.put(n);
			index.put(n,ind);
			ind++;
		}
	}

	public static int getNumberOfProcessors(){return MAX_THREADS;}

	public  long getResults(){
		return results.get();
	}

	public HashMap<Integer,Integer> getIndexMap(){
		return index;
	}

	public Double call(){
		double fail=0.0;
		Integer taskInd;
		//long id=Thread.currentThread().getId();
		while ((taskInd = tasks.poll()) != null ){ 
			//System.out.println("Thread [" + id + "]: Processing: " + taskInd);
			//results[index.get(taskInd)]=compute(taskInd).clone();
			fail = compute(taskInd);
		}
		return Double.valueOf(fail);
	}


	public double compute(Integer record){
		//Execute command
		double y = (Double) null; //standard error
			y=0;
		return y;
	}

	public long start() throws Exception{
		
		DistanceMatrix map=null;
		int threads= MAX_THREADS;
		threads=(threads > MAX_THREADS) ? MAX_THREADS:threads;
		threads=(threads < 1) ? MAX_THREADS:threads;
		
		ArrayList<Future<Double>> futures=
				new ArrayList<Future<Double>>(threads);
		ExecutorService executor= 
				Executors.newFixedThreadPool(threads);

		double fail=0;
		try {
			map = new DistanceMatrix(database);

			for(int i=0;i<threads;i++){
				Future<Double> future= executor.submit(map);
				futures.add(future);
			}
			for( Future<Double> future: futures){
				fail =future.get();
				if(fail>0)
					System.err.println("Future computation failed:  " + future.toString());
			}
			results.addAndGet((long) fail);
		} catch (InterruptedException e1) {
			e1.printStackTrace();
		} catch (ExecutionException e) {
			e.printStackTrace();
		} 
		executor.shutdown();
		//System.out.println("!!Done: Processed records: " + results.length);
		return results.get();
	}

}
