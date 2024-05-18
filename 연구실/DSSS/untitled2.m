% 데이터 생성
x = linspace(0, 2*pi, 100);
y1 = sin(x);
y2 = cos(x);

% 그래프 그리기
figure;
plot(x, y1, 'r', 'LineWidth', 2); % 빨간색 선
hold on;
plot(x, y2, 'b', 'LineWidth', 2); % 파란색 선
hold off;

% 그래프 설정 (배경색을 흰색으로 설정)
set(gcf, 'Color', 'w'); % figure 배경을 흰색으로 설정
set(gca, 'Color', 'w'); % axes 배경을 흰색으로 설정

% 그래프 저장
print('clean_graph', '-dpng', '-r300');
