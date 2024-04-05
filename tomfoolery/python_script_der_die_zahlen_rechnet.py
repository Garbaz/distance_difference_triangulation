PERIOD = 1.0
DELAY = 0.05678

def send(prev_send, prev_rec):
    if prev_rec > prev_send:
        return (prev_send + prev_rec)/2 + 1.5 * PERIOD
    else:
        print("!!!")
        return prev_send + 2 * PERIOD
    
    # return prev_send + 2 * PERIOD


send_a = 0.0
rec_a = 0.0
send_b = 0.1234
rec_b = 0.0

while True:
    send_a = send(send_a, rec_a)
    rec_b = send_a + DELAY
    send_b = send(send_b, rec_b)
    rec_a = send_b + DELAY
    print(f"A: {send_a:<16} | {rec_a:<16}")
    print(f"B: {send_b:<16} | {rec_b:<16}")
    print(f"{send_b - send_a}")
    input()
