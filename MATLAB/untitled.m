clear all
close all

% Video dosyasını yükle
video_file = 'vidyo.mp4'; % Video dosya adını buraya yazın
v = VideoReader(video_file);

% Çıktı video ayarları
output_video = VideoWriter('crack_analysis_output.mp4', 'MPEG-4');
output_video.FrameRate = v.FrameRate;
open(output_video);

% Sonuçları saklamak için değişkenler
frame_count = 0;
total_cracks = [];
avg_distances = [];
min_distances = [];
max_distances = [];
crack_areas = [];
crack_ratios = []; % Çatlak/yüzey oranları

% Filtreleme parametreleri
G = fspecial('gaussian',[9 9],3);
Sh = 9;  % Yüksek eşik
Sb = 6;  % Düşük eşik
se = strel('square',4);

fprintf('Video işleme başlıyor...\n');

% Video frame'lerini işle
while hasFrame(v)
    frame_count = frame_count + 1;
    I = readFrame(v);
    
    % Orijinal kod mantığını frame'e uygula
    Indg = rgb2gray(I);
    [l,c] = size(Indg);
    
    % Gaussian filtreleme
    IG = imfilter(Indg, G, 'symmetric', 'same', 'conv');
    
    % Gradyan hesaplama
    IX = zeros(l,c);
    IY = zeros(l,c);
    IXY = zeros(l,c);
    
    for i = 2:l-1
        for j = 2:c-1
            IX(i,j) = -(-double(IG(i-1,j-1))-2*double(IG(i,j-1))-double(IG(i+1,j-1)) + ...
                double(IG(i-1,j+1))+2*double(IG(i,j+1))+double(IG(i+1,j+1)))/8;
            IY(i,j) = -(-2*double(IG(i-1,j))+2*double(IG(i+1,j))-double(IG(i-1,j-1))+ ...
                double(IG(i+1,j-1))-double(IG(i-1,j+1))+double(IG(i+1,j+1)))/8;
            IXY(i,j) = sqrt(IX(i,j)^2 + IY(i,j)^2);
        end
    end
    
    % Binary görüntü oluşturma
    Ibin = zeros(l,c);
    for i = 1:l
        for j = 1:c
            if IXY(i,j) >= Sh
                Ibin(i,j) = 255;
            elseif IXY(i,j) < Sb
                Ibin(i,j) = 0;
            end
        end
    end
    
    % Bağlantı işlemi
    for p = 1:20
        for i = 2:l-1
            for j = 2:c-1
                if IXY(i,j) < Sh && IXY(i,j) >= Sb
                    vect = [Ibin(i-1,j-1:j+1) Ibin(i,j-1) Ibin(i,j+1) Ibin(i+1,j-1:j+1)];
                    if max(vect) == 255
                        Ibin(i,j) = 255;
                    else
                        Ibin(i,j) = 0;
                    end
                end
            end
        end
    end
    
    % Bölge analizi ve gürültü temizleme
    [L, num] = bwlabel(Ibin, 8);
    t = regionprops(L, 'area');
    taille = zeros(num,1);
    for g = 1:num
        taille(g,1) = t(g).Area;
    end
    
    Ifer = zeros(l,c);
    for i = 1:l
        for j = 1:c
            if Ibin(i,j) == 0
                Ifer(i,j) = 0;
            else
                for w = 1:num
                    if L(i,j) == w
                        if taille(w) <= 500
                            Ifer(i,j) = 0;
                        else
                            Ifer(i,j) = Ibin(i,j);
                        end
                    end
                end
            end
        end
    end
    
    % Morfolojik closing işlemi
    Ifer2 = imclose(Ifer, se);
    
    % Final bölge analizi
    [L2, num2] = bwlabel(Ifer2, 8);
    t2 = regionprops(L2, 'area');
    taille2 = zeros(num2,1);
    for g = 1:num2
        taille2(g,1) = t2(g).Area;
    end
    
    % Sonuç görüntüsü oluşturma (çatlakları kırmızı ile işaretle)
    Ifin = I;
    for i = 1:l
        for j = 1:c
            if Ifer2(i,j) ~= 0
                Ifin(i,j,1) = 255; % Kırmızı kanal
                Ifin(i,j,2) = 0;   % Yeşil kanal
                Ifin(i,j,3) = 0;   % Mavi kanal
            end
        end
    end
    
    % Çatlak merkezlerini hesapla
    stats = regionprops(L2, 'Centroid');
    centroids = cat(1, stats.Centroid);
    
    % Merkezleri görüntüye ekle
    if ~isempty(centroids)
        for k = 1:size(centroids,1)
            x = round(centroids(k,1));
            y = round(centroids(k,2));
            if x >= 3 && x <= c-2 && y >= 3 && y <= l-2
                % Sarı çarpı işareti ekle
                Ifin(y-2:y+2, x, :) = repmat([255 255 0], 5, 1); % Dikey çizgi
                Ifin(y, x-2:x+2, :) = repmat(reshape([255 255 0], 1, 1, 3), 1, 5, 1); % Yatay çizgi
            end
        end
    end
    
    % Frame bilgilerini kaydet
    total_cracks = [total_cracks; num2];
    total_area = sum(taille2);
    crack_areas = [crack_areas; total_area];
    
    % Çatlak/yüzey oranını hesapla
    total_surface_area = l * c; % Toplam piksel sayısı
    crack_ratio = (total_area / total_surface_area) * 100; % Yüzde olarak
    crack_ratios = [crack_ratios; crack_ratio];
    
    % Çatlaklar arası mesafe analizi
    if size(centroids,1) > 1
        distances = pdist(centroids);
        avg_distances = [avg_distances; mean(distances)];
        min_distances = [min_distances; min(distances)];
        max_distances = [max_distances; max(distances)];
    else
        avg_distances = [avg_distances; NaN];
        min_distances = [min_distances; NaN];
        max_distances = [max_distances; NaN];
    end
    
    % Frame üzerine bilgi metni ekle
    text_info = sprintf('Frame: %d | Çatlaklar: %d | Alan: %.0f px² | Oran: %.3f%%', ...
        frame_count, num2, total_area, crack_ratio);
    text_img = insertText(Ifin, [10 10], text_info, 'FontSize', 12, 'BoxColor', 'white', 'TextColor', 'black');
    
    % İşlenmiş frame'i çıktı videosuna yaz
    writeVideo(output_video, text_img);
    
    % İlerleme göstergesi
    if mod(frame_count, 30) == 0
        fprintf('İşlenen frame sayısı: %d\n', frame_count);
    end
end

% Video dosyalarını kapat
close(output_video);
fprintf('Video işleme tamamlandı. Toplam frame: %d\n', frame_count);

% Sonuçları analiz et ve görselleştir
fprintf('\n=== VİDEO ÇATLAK ANALİZİ SONUÇLARI ===\n');
fprintf('Toplam frame sayısı: %d\n', frame_count);
fprintf('Ortalama çatlak sayısı: %.2f\n', mean(total_cracks));
fprintf('Maksimum çatlak sayısı: %d\n', max(total_cracks));
fprintf('Minimum çatlak sayısı: %d\n', min(total_cracks));
fprintf('Ortalama çatlak alanı: %.2f px²\n', mean(crack_areas));
fprintf('Ortalama çatlak/yüzey oranı: %.4f%% (%.6f)\n', mean(crack_ratios), mean(crack_ratios)/100);
fprintf('Maksimum çatlak/yüzey oranı: %.4f%%\n', max(crack_ratios));
fprintf('Minimum çatlak/yüzey oranı: %.4f%%\n', min(crack_ratios));

% Çatlak şiddeti kategorileri
severe_frames = sum(crack_ratios > 1.0); % %1'den fazla çatlak
moderate_frames = sum(crack_ratios > 0.5 & crack_ratios <= 1.0); % %0.5-1 arası
mild_frames = sum(crack_ratios > 0.1 & crack_ratios <= 0.5); % %0.1-0.5 arası
minimal_frames = sum(crack_ratios <= 0.1); % %0.1'den az

fprintf('\n=== ÇATLAK ŞİDDETİ KATEGORİLERİ ===\n');
fprintf('Ciddi çatlak (>%%1.0): %d frame (%.1f%%)\n', severe_frames, (severe_frames/frame_count)*100);
fprintf('Orta seviye (%%0.5-1.0): %d frame (%.1f%%)\n', moderate_frames, (moderate_frames/frame_count)*100);
fprintf('Hafif çatlak (%%0.1-0.5): %d frame (%.1f%%)\n', mild_frames, (mild_frames/frame_count)*100);
fprintf('Minimal çatlak (≤%%0.1): %d frame (%.1f%%)\n', minimal_frames, (minimal_frames/frame_count)*100);

% NaN değerleri olmayan mesafeleri filtrele
valid_avg_dist = avg_distances(~isnan(avg_distances));
valid_min_dist = min_distances(~isnan(min_distances));
valid_max_dist = max_distances(~isnan(max_distances));

if ~isempty(valid_avg_dist)
    fprintf('Ortalama çatlaklar arası mesafe: %.2f px\n', mean(valid_avg_dist));
    fprintf('Minimum çatlaklar arası mesafe: %.2f px\n', mean(valid_min_dist));
    fprintf('Maksimum çatlaklar arası mesafe: %.2f px\n', mean(valid_max_dist));
end

% Grafikleri oluştur
figure;
subplot(2,2,1);
plot(1:frame_count, total_cracks, 'b-', 'LineWidth', 2);
title('Frame Başına Çatlak Sayısı');
xlabel('Frame Numarası');
ylabel('Çatlak Sayısı');
grid on;

subplot(2,2,2);
plot(1:frame_count, crack_areas, 'r-', 'LineWidth', 2);
title('Frame Başına Toplam Çatlak Alanı');
xlabel('Frame Numarası');
ylabel('Alan (px²)');
grid on;

subplot(2,2,3);
if ~isempty(valid_avg_dist)
    plot(1:length(valid_avg_dist), valid_avg_dist, 'g-', 'LineWidth', 2);
    title('Çatlaklar Arası Ortalama Mesafe');
    xlabel('Frame Numarası (Geçerli)');
    ylabel('Mesafe (px)');
    grid on;
end

subplot(2,2,4);
histogram(total_cracks, 'BinWidth', 1, 'FaceColor', 'cyan');
title('Çatlak Sayısı Dağılımı');
xlabel('Çatlak Sayısı');
ylabel('Frame Sayısı');
grid on;

% Sonuçları Excel dosyasına kaydet (opsiyonel)
results_table = table((1:frame_count)', total_cracks, crack_areas, avg_distances, min_distances, max_distances, ...
    'VariableNames', {'Frame', 'CrackCount', 'TotalArea', 'AvgDistance', 'MinDistance', 'MaxDistance'});

try
    writetable(results_table, 'crack_analysis_results.xlsx');
    fprintf('\nSonuçlar crack_analysis_results.xlsx dosyasına kaydedildi.\n');
catch
    fprintf('\nExcel dosyası kaydedilemedi. Sadece workspace''te results_table değişkeni mevcut.\n');
end

fprintf('\nAnaliz tamamlandı!\n');
fprintf('Çıktı videosu: crack_analysis_output.mp4\n');