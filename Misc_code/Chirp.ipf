#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later



Function ZAP_analysis(ctrlName) : ButtonControl
	String  ctrlName
	
	execute/z "SetDataFolderFunction()"
	NVAR num_channel=root:num_channel///
	string/g tracelist, tracename
	newdatafolder/o :ZAPanalysis 
	tracelist=tracenamelist("",",",1)
	tracename=Stringfromlist(0,tracelist,",")
	duplicate/O $tracename, :ZAPanalysis:tmpwave
	
	
	setdatafolder :ZAPanalysis
	variable z=num_channel-1
	variable i=num_channel-1//counter for current injection
	variable j=num_channel   //counter for voltage response
	variable k=0 //counter for number of channels
	variable sampling_rate
	make/o/n=(num_channel) MaxZ, Fr, cutoff, zap_tau, phase0, Iphase,zap_tau,RQ 
								
	wave tmpwave //duplicate each channel 
	
	do
		duplicate/o/R=[,][i-z,i-z] tmpwave, current ///This is the bulk of the ZAP and ZPP code starting here
		duplicate/o/R=[,][j,j]  tmpwave, Vm 
		duplicate/o Vm, Vmx
		Vmx=x
		sampling_rate=Vmx[1]-Vmx[0]
		sampling_rate=1/sampling_rate
		sampling_rate*=1000
		
		print "sampling rate is",sampling_rate
	
	FFT/DEST=current_FFT current									//FFT of the current and response
	FFT/DEST=Vm_FFT Vm
	
	Duplicate/o/c Vm_FFT, ZAP
	ZAP=(Vm_FFT)/(current_FFT)
	ZAP*=1000
														
		wavestats/q/c=1 Vm_FFT
		make/o/n=(V_npnts) ZapOutput, Zapimg, ZAPreal, tmp_phase, phase
		ZapOutput=sqrt(magsqr(ZAP))
		ZAPimg=imag(ZAP)
		ZAPreal=real(ZAP)
		tmp_phase=ZAPimg/ZAPreal
		phase=atan(tmp_phase)*-1
		SetScale/I x 0, sampling_rate/2,"", ZapOutput//adjust this to match sampling rate
		SetScale/I x 0, sampling_rate/2,"", phase///
		
			
	duplicate/o ZAPOutput, ZAP_smooth
	smooth/m=0 20, ZAP_smooth								//smooths the ZAP
	duplicate/o phase, phase_smooth
	smooth/m=0 20, phase_smooth
	
	wavestats/q/r=(.5, 40) ZAP_smooth
	Fr[k]=v_maxloc
	MaxZ[k]=V_max
	RQ[k]=MaxZ[k]
	RQ[k]/=zap_smooth(1)
	zap_tau[k]=MaxZ[k]
	zap_tau[k]/=zap_smooth(1)
	duplicate/o ZAP_smooth, ZAP_norm								//normalizes the ZAP
	ZAP_norm/=MaxZ[k]
		findlevel/EDGE=2/q/b=5/q/r=(.5,40), ZAP_norm, .5^(.5)
		cutoff[k]=v_levelx

	zap_tau[k]=(1/(cutoff[k]*(2*pi)))*1000 // estimate time constant of underlying filter
	
	
	phase_smooth*=(180/pi) //convert phase to degrees
	findlevel/q/r=(.5,40) phase_smooth, 0
	phase0[k]=v_levelx
	If  (waveexists(phase0))
	
		Iphase[k]=area(phase_smooth, .5, phase0[0])
		else 
		endif
	
	
	Display ZAP_smooth
	SetAxis left 0,300
	SetAxis bottom .5, 40
	ModifyGraph log(bottom)=1
	
	Display ZAP_norm
	SetAxis left 0, 1.25
	SetAxis bottom, .5, 40
	ModifyGraph log(bottom)=1
	
	Display phase_smooth
	SetAxis bottom .5, 40
	ModifyGraph log(bottom)=1
	
		i+=1
		j+=1
		k+=1
		
	
	while (k<num_channel)	
	
	make/o/n=1 Zap_first,Zap_last
	NVAR first=root:first_trace
 	Zap_first=first
 	NVAR last=root:last_trace
 	Zap_last=last
	edit root:cell_name, Fr, RQ , cutoff, zap_tau,Iphase, phase0, MaxZ,Zap_first,Zap_last
	
 end	