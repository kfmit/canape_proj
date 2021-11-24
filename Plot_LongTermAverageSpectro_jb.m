
function Plot_LongTermAverageSpectro_jb(A,FrequencyVector,TimeVector,FeatureMatrix,Opt, db_lim_plot)

UseTimeStamp= 1;

% Long-term spectrogram

clf; 
imagesc(TimeVector,FrequencyVector,A.');
colormap(1-gray)
set(gca,'YDir','normal')
ylabel('Frequency [ Hz ]')
ylabel(colorbar,'PSD [ dB re 1 \muPa^2 Hz^-^1 ]','fontname','arial','fontsize',14)

if ~isempty(db_lim_plot)
    caxis(db_lim_plot)
end

[LineSpecificationMatrix,~,ColorMatrix] = LineAndColorSpecificationMatrix(size(FeatureMatrix,2)); 

hold on,
for cc= 1:size(FeatureMatrix,2)  
        
    Val_Desc = FeatureMatrix(:,cc);
 
    if length(find(Val_Desc<0))>0
        Val_Desc = Val_Desc - nanmin(Val_Desc);        
    end
    
    if Opt.MV_Apply_MedFilt>0
        Val_Desc = medfilt1(Val_Desc,Opt.MV_Apply_MedFilt,[],1,'omitnan');
    end
    
%     plot(TimeVector,...
%         Opt.OffsetFreqDescriptors + FrequencyVector(end)/3 * abs(Val_Desc) / max(abs(Val_Desc)),...
%         '-','linewidth',1.2,'color',ColorMatrix(cc,:),'markersize',8);
    plot(TimeVector,...
        Opt.OffsetFreqDescriptors + FrequencyVector(end)/3 * abs(Val_Desc) / max(abs(Val_Desc)),...
        LineSpecificationMatrix{cc},'linewidth',1.2,'color',ColorMatrix(cc,:),'markersize',8);

    
    Opt.OffsetFreqDescriptors = Opt.OffsetFreqDescriptors + Opt.InterFreqDescriptors;

end
legend(Opt.AuxVariableNames)
    
if UseTimeStamp   
   PutTimeStamp(TimeVector,Opt) 
end

end

