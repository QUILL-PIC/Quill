#include <iostream>
#include <string.h>
#include <sys/times.h>
#include <unistd.h>
#include <algorithm>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>

#include "main.h"

extern ddi* p_last_ddi;

void* listen_for_param_updates(void*)
{
    cout << "[listener] thread started" << endl;
    
    int listenfd = 0, connfd = 0;
    int n_bytes = 0;
    char recvBuff[1024];
    struct sockaddr_in serv_addr;
    int server_port;
    struct sockaddr_storage client_addr;
    int client_port;
    socklen_t len;
    char ipstr[INET6_ADDRSTRLEN];

    server_port = 5700 + int(200.0 * rand() / RAND_MAX);

    listenfd = socket(AF_INET, SOCK_STREAM, 0);
    memset(&serv_addr, '0', sizeof(serv_addr));
    memset(recvBuff, '0' ,sizeof(recvBuff));

    serv_addr.sin_family = AF_INET;
    serv_addr.sin_addr.s_addr = htonl(INADDR_ANY);
    serv_addr.sin_port = htons(server_port); 

    bind(listenfd, (struct sockaddr*)&serv_addr, sizeof(serv_addr)); 

    listen(listenfd, 10); 
    
    cout << "[listener] Started listening for param updates. Port = " << server_port << endl;
    cout << "[listener] Connect clients using nc: nc [hostname] " << server_port << endl;

    while(1)
    {
        len = sizeof client_addr;
        connfd = accept(listenfd, (struct sockaddr*)&client_addr, &len); 
        if (client_addr.ss_family == AF_INET) 
        {
            struct sockaddr_in *s = (struct sockaddr_in *)&client_addr;
            client_port = ntohs(s->sin_port);
            inet_ntop(AF_INET, &s->sin_addr, ipstr, sizeof ipstr);
        }
        cout << "[listener] Received new connection from " << ipstr << " port = " << client_port << endl;

        while((n_bytes = read(connfd, recvBuff, sizeof(recvBuff)-1)) > 0)
        {
            recvBuff[n_bytes] = 0;
 
            if(n_bytes < 0)
            {
                cout << "[listener] Error - unable to read bytes from client" << endl;
            }
            else
            {
                cout << "[listener] " << n_bytes << " bytes received: " << recvBuff << endl;
            }
            
            string received(recvBuff);
            //cout << received.length() << "String = " << received << " Last byte: " << int(received[3]) << endl;
            
            // remove whitespaces
            received.erase(remove(received.begin(), received.end(), ' '), received.end()); 
            int pos = received.find("t_end=");
            double t_end_received = 0.0;
            if (pos != string::npos)
            {
                const char* a = received.substr(6).c_str();
                //cout << "to parse: " << a << endl;
                t_end_received = atof(a);
            }
            
            if (t_end_received != 0.0)
            {
                cout << "[listener] Parsed t_end: " << t_end_received << endl;
                if (t_end_received > 0.0)
                {
                    cout << "[listener] Updating p_last_ddi->t_end: [" << p_last_ddi->t_end/2/PI << "] => [" << t_end_received << "]" << endl;
                    p_last_ddi->t_end = t_end_received * 2 * PI; // need synchronization?
                }
            }
        }

        close(connfd);
        sleep(1);
    }
}
