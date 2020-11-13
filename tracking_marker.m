clc;close all;imtool close all;clear;

% トラッキングする最大数
numTracks = 5;

% 読み込むビデオのフレームの最大数
numVideoFrames = 200;

% Stop ボタン表示
runLoop=true;
sz = get(0,'ScreenSize');
figure('MenuBar','none','Toolbar','none','Position',[20 sz(4)-100 100 70])
uicontrol('Style', 'pushbutton', 'String', 'Stop',...
    'Position', [20 20 80 40], 'Callback', 'runLoop=false;');

% 動画読込みのオブジェクト生成とフレームレート取得
videoReader = vision.VideoFileReader("D:\pictures\DSCF4161_Trim_Trim_strict.mp4");

% サンプルレートを取得
S = videoReader.info;
frameRate = S.VideoFrameRate;

% ビデオサイズを取得
videoSize = S.VideoSize;

% 動画書き出しのオブジェクト生成
videoWriter = vision.VideoFileWriter('result.mp4','FileFormat','MPEG4','FrameRate',frameRate);

% 動画表示用オブジェクトの生成
videoPlayer = vision.VideoPlayer('Position', [ 20, sz(4)-videoSize(2)-250, videoSize+50]);

% ブロブ解析用のオブジェクトの作成（中心点・境界ボックス）
blobAnalyser = vision.BlobAnalysis( 'AreaOutputPort',false, ...
    'CentroidOutputPort', true,  'BoundingBoxOutputPort', true, ...
    'MinimumBlobArea', 100, 'ExcludeBorderBlobs',true);

% トラック（各トラックに、各動いている物体の情報を格納）の初期化
tracks = struct(...
    'id',           {}, ...       % ID番号
    'bbox',         {}, ...       %現フレーム内の境界ボックス（表示用）
    'kalmanFilter', {}, ...       %この物体トラッキング用のカルマンフィルター
    'age',               {}, ...  %最初に検出されてからのフレーム数
    'totalVisibleCount', {}, ...      %検出された全フレーム数           => これが閾値を超えたら表示する
    'consecutiveInvisibleCount', {}); %連続した非検出フレーム数         => これが閾値を超えたらこのトラックを削除

nextId = 1; % ID of the next track
centroidLog = nan(numTracks,2*numVideoFrames);
frameCount = 1;

% 1フレーム目を読み込み
frame = step(videoReader);

% 色のしきい値アプリを使ってマスク関数を作成
colorThresholder(frame);
% 色空間はHSVを選択
% 緑色のマーカーだけが検出されるようにしきい値を決める
% しきい値が決まったら「エクスポート」⇒「関数のエクスポート」で
% createMask.m関数を生成し、保存する

% ビデオの読み込みをリセット
videoReader.reset;

while (runLoop && ~isDone(videoReader)) %runLoopはフラグ、isDone()はvideoReaderが最後までいったらtrue
    frame = step(videoReader);     % 1フレーム読み込み

    % フレーム内の物体を検出：動いている物体(領域)の検出：maskは1が前景・0が背景
    % 下記のcreateMask関数はcolorThresholderアプリケーションを使って作成
    mask = createMask(frame);

    % 中心点・境界ボックスの検出
    [centroids, bboxes] = step(blobAnalyser, mask);

    % 現フレーム内での位置を前フレームから予測
    % カルマンフィルタを用いて、前フレームまでに検出済みの各物体(トラック)の位置予測
    for i = 1:length(tracks)
        bbox = tracks(i).bbox;    % 前フレームでの境界ボックス
        % 現フレームでの位置の予測
        predictedPosition = int32(predict(tracks(i).kalmanFilter));
        % 境界ボックスの中心を、予測重心位置へ調整
        tracks(i).bbox = [predictedPosition - bbox(3:4)/2, bbox(3:4)];

    end

    % 既検出物体(トラック)に、検出された物体を対応付け
    % コストの計算：既検出物体の予測位置と各検出物体の距離
    cost = zeros(length(tracks), size(centroids, 1));         % トラック数 x 検出数
    for i = 1:length(tracks)
        cost(i, :) = distance(tracks(i).kalmanFilter, centroids);
    end

    % 各トラックに検出された物体を割当 (ハンガリアンアルゴリズム：コストの合計が最小になるように）
    % assignment:トラック番号と検出番号の組が代入される
    %costOfNonAssignment = 20;     %小さいと既存トラックにアサインされないもの増える=>新しいトラックが多く生成される
    costOfNonAssignment = 50;     %小さいと既存トラックにアサインされないもの増える=>新しいトラックが多く生成される
    [assignments, unassignedTracks, unassignedDetections] = ...
        assignDetectionsToTracks(cost, costOfNonAssignment);

    % 対応する物体が見つかったトラックの情報更新
    for i = 1:size(assignments, 1)
        trackIdx = assignments(i, 1);           %トラックの番号
        detectionIdx = assignments(i, 2);       %検出物体の番号
        centroid = centroids(detectionIdx, :);  %その検出された物体の中心点を抽出（blobAnalyser.step(mask)の結果を使用）
        bbox = bboxes(detectionIdx, :);         %その検出された物体の境界ボックスを抽出（blobAnalyser.step(mask)の結果を使用)
        % 検出された位置と予測値を用いて、現在位置を推定
        correct(tracks(trackIdx).kalmanFilter, centroid);

        % 各トラック情報を更新：境界ボックス、age、totalVisibleCount、consecutiveInvisibleCount
        tracks(trackIdx).bbox = bbox;
        tracks(trackIdx).age = tracks(trackIdx).age + 1;
        tracks(trackIdx).totalVisibleCount = ...
            tracks(trackIdx).totalVisibleCount + 1;
        tracks(trackIdx).consecutiveInvisibleCount = 0;
    end

    % 対応する物体が見つからなかったトラックの情報更新
    %           age、consecutiveInvisibleCount
    for i = 1:length(unassignedTracks)
        ind = unassignedTracks(i);
        tracks(ind).age = tracks(ind).age + 1;
        tracks(ind).consecutiveInvisibleCount = ...
            tracks(ind).consecutiveInvisibleCount + 1;
    end

    % 見失ったトラックを消去
    % 連続5フレーム以上不観測でトラックから消去
    if ~isempty(tracks)
        % compute the fraction of the track's age for which it was visible
        ages = [tracks(:).age];
        totalVisibleCounts = [tracks(:).totalVisibleCount];
        visibility = totalVisibleCounts ./ ages;     % 全age中で観測されたフレームの割合

        % find the indices of 'lost' tracks
        lostInds = (ages < 4 & visibility < 0.6) | ... 
            ([tracks(:).consecutiveInvisibleCount] >= 40);   %もとは5, 30フレーム以上連続不観測で消去

        tracks = tracks(~lostInds);   %見失ったトラックの消去
    end

    % 新たに見つかった物体に対し、新しいトラックを生成
    %      (ここではアサインされなかったトラックを、新しいトラックとする)
    centroids1 = centroids(unassignedDetections, :);           % Nx2 double
    bboxes = bboxes(unassignedDetections, :);                 % Nx4 int32

    for i = 1:size(centroids1, 1)

        centroid = centroids1(i,:);
        bbox = bboxes(i, :);

        % あたらしい物体一つに対して、カルマンフィルターを1つ生成
        % http://www.mathworks.com/videos/introduction-to-kalman-filters-for-object-tracking-79674.html
        kalmanFilter = configureKalmanFilter('ConstantVelocity', ...
            centroid, [2, 5], [1, 2.5], 10); % [200, 50], [100, 25], 100);

        % 新しいトラックの生成
        newTrack = struct(...
            'id', nextId, ...
            'bbox', bbox, ...
            'kalmanFilter', kalmanFilter, ...
            'age', 1, ...
            'totalVisibleCount', 1, ...
            'consecutiveInvisibleCount', 0);

        % 生成した新しいトラックを、トラックの配列の最後に追加
        tracks(end + 1) = newTrack;

        % nextIdを1つ増やす
        nextId = nextId + 1;
    end


    % 結果の表示
    frame = im2uint8(frame);

    if ~isempty(tracks)

        reliableTrackInds = ...
            [tracks(:).totalVisibleCount] > 0;
        reliableTracks = tracks(reliableTrackInds);

        if ~isempty(reliableTracks)
            % 境界の四角枠の座標の抽出
            bboxes = cat(1, reliableTracks.bbox);

            % idsの取得
            ids = int32([reliableTracks(:).id]);

            % create labels for objects indicating the ones for
            % which we display the predicted rather than the actual
            % location
            labels = cellstr(int2str(ids'));
            predictedTrackInds = ...
                [reliableTracks(:).consecutiveInvisibleCount] > 0;
            isPredicted = cell(size(labels));
            isPredicted(predictedTrackInds) = {' predicted'};
            labels = strcat(labels, isPredicted);

            % RGBのフレームの中に、四角枠を描画
            frame = insertObjectAnnotation(frame, 'rectangle', ...
                bboxes, labels);
        end
    end

    % 奇跡の保存と可視化
    ids = int32([tracks(:).id]);
    bboxTracked = cat(1,tracks.bbox);
    centroidLog(ids,(frameCount-1)*2+1:(frameCount-1)*2+2) = bboxTracked(:,1:2);

    % 表示
    videoPlayer.step(frame);

    % 動画に書き出し
    step(videoWriter,frame);

    % フレームのカウントを1増やす
    frameCount = frameCount + 1;

    % プッシュボタンのイベントの確認
    drawnow;
end
release(videoReader);
release(videoPlayer);
release(videoWriter);

%% 軌跡を可視化
figure;
image(frame,'AlphaData',0.3);
axis image;
hold on;
x = centroidLog(:,1:2:end)';
y = centroidLog(:,2:2:end)';
plot(x,y);
title('軌跡'); xlabel('X軸位置(pixel)'); ylabel('Y軸位置(pixel)');
axis ij;

%% 位置の時間変化の可視化
figure;
t = ((0:(size(centroidLog,2)/2-1))/frameRate)';
subplot(2,1,1),plot(t,x);
title('X軸'); xlabel('時間(sec)'); ylabel('X軸位置(pixel)');
subplot(2,1,2),plot(t,y);
title('Y軸'); xlabel('時間(sec)'); ylabel('Y軸位置(pixel)');

%% CSVファイルとして書き出し
% table型に変換
tbl = table(t,x,y);
disp(tbl(1:10,:));

% CSVとして書き出し
writetable(tbl,'result.csv');

% CSVを読み出して確認
tbl2 = readtable('result.csv');
disp(tbl2(1:5,[1:4 102:104]));