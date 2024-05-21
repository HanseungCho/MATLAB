clc; clear; close all
N = 100;% Number of users
G = linspace(0, 8, 100);

throughput_aloha = G .* exp(-2*G);
throughput_slotted = G .* exp(- G);

sim_throughput_aloha = zeros(size(G));
sim_throughput_slotted = zeros(size(G));
iterations = 10000;
%상황: 첫 packet이 전송 된 상황에서 그 다음에 발생될 packet과의 collision simulation
for i = 1:length(G)
    p = G(i)/N;  % 한 유저당 packet보낼 확률

    % ALOHA Simulation
    success_aloha = 0;
    for j = 1:iterations
        for k=1:20 %2T(20 mini slots)
            sent_packets0(k) = binornd(N, p/10);% slot time T를 10개의 mini slot으로 모델링
        end %2T에서 각 slot에서 100명의 유저가 보낼 packet 수 생성
        if sum(sent_packets0) == 1 %2T 안에서 packet이 하나만 보내진 사건을 카운트
            success_aloha = success_aloha + 1;
        end
    end 
    sim_throughput_aloha(i) = success_aloha / iterations;
    %T동안의 Throughput을 구해야하므로 1/2적용
    sim_throughput_aloha(i)=sim_throughput_aloha(i)/2;

    % Slotted ALOHA Simulation
    %한 슬롯에 packet이 하나만 전송되는 확률
    success_slotted = 0;
    for j = 1:iterations
        sent_packets1 = binornd(N, p);
        if sent_packets1 == 1 
            success_slotted = success_slotted + 1;
        end
    end
    sim_throughput_slotted(i) = success_slotted / iterations;
end

% Plot results
figure;
hold on;
plot(G, throughput_aloha);
plot(G, throughput_slotted);
plot(G, sim_throughput_aloha, 'b--', 'LineWidth', 1.5);
plot(G, sim_throughput_slotted, 'k--', 'LineWidth', 1.5);
xlabel('Traffic load G');
ylabel('Throughput');
legend('Theoretical ALOHA', 'Theoretical Slotted ALOHA', 'Simulated ALOHA', 'Simulated Slotted ALOHA');
grid on;
hold off;
